import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
from collections import defaultdict
import pandas as pd
from matplotlib.figure import Figure

class PCRGelApp:
    """Aplicación interactiva para simulación de PCR y electroforesis en gel"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("In Silico PCR Gel Simulator")
        self.root.geometry("1200x700")
        
        # Datos de la aplicación
        self.samples = []  # Lista de (nombre, tamaño)
        
        # Escaleras disponibles
        self.ladders = {
            '1 kb': [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 750, 500, 250],
            '100 bp': [1500, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
            'Lambda': [23130, 9416, 6557, 4361, 2322, 2027, 564, 125],
            'Low Range': [1000, 800, 600, 500, 400, 300, 200, 100, 80, 60, 40, 20]
        }
        
        # Configuración del gel
        self.gel_percentage = tk.DoubleVar(value=1.0)
        self.voltage = tk.IntVar(value=100)
        self.runtime = tk.IntVar(value=60)
        self.selected_ladder = tk.StringVar(value='1 kb')
        
        # Crear la interfaz
        self.create_widgets()
        
        # Inicializar el gráfico
        self.update_gel()
    
    def create_widgets(self):
        """Crea todos los widgets de la interfaz"""
        
        # Frame principal dividido en paneles
        main_panel = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_panel.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Panel izquierdo - Controles
        left_panel = ttk.Frame(main_panel)
        main_panel.add(left_panel, weight=1)
        
        # Panel derecho - Gráfico
        right_panel = ttk.Frame(main_panel)
        main_panel.add(right_panel, weight=2)
        
        # ===== PANEL IZQUIERDO =====
        # Título
        title_label = ttk.Label(left_panel, text="Gel Configuration", 
                                font=('Arial', 14, 'bold'))
        title_label.pack(pady=10)
        
        # Separador
        ttk.Separator(left_panel, orient='horizontal').pack(fill='x', padx=10, pady=5)
        
        # Configuración del gel
        config_frame = ttk.LabelFrame(left_panel, text="Gel Parameters", padding=10)
        config_frame.pack(fill='x', padx=10, pady=5)
        
        # Porcentaje de agarosa
        ttk.Label(config_frame, text="Agarose (%):").grid(row=0, column=0, sticky='w', pady=2)
        agarose_combo = ttk.Combobox(config_frame, textvariable=self.gel_percentage, 
                                     values=[0.8, 1.0, 1.2, 1.5, 2.0], width=10)
        agarose_combo.grid(row=0, column=1, pady=2, padx=5)
        agarose_combo.bind('<<ComboboxSelected>>', lambda e: self.update_gel())
        
        # Voltaje
        ttk.Label(config_frame, text="Voltage (V):").grid(row=1, column=0, sticky='w', pady=2)
        voltage_spin = ttk.Spinbox(config_frame, from_=50, to=200, increment=10, 
                                   textvariable=self.voltage, width=10)
        voltage_spin.grid(row=1, column=1, pady=2, padx=5)
        voltage_spin.bind('<KeyRelease>', lambda e: self.update_gel())
        
        # Tiempo
        ttk.Label(config_frame, text="Run Time (min):").grid(row=2, column=0, sticky='w', pady=2)
        runtime_spin = ttk.Spinbox(config_frame, from_=20, to=120, increment=5, 
                                   textvariable=self.runtime, width=10)
        runtime_spin.grid(row=2, column=1, pady=2, padx=5)
        runtime_spin.bind('<KeyRelease>', lambda e: self.update_gel())
        
        # Ladder
        ttk.Label(config_frame, text="DNA Ladder:").grid(row=3, column=0, sticky='w', pady=2)
        ladder_combo = ttk.Combobox(config_frame, textvariable=self.selected_ladder,
                                     values=list(self.ladders.keys()), width=10)
        ladder_combo.grid(row=3, column=1, pady=2, padx=5)
        ladder_combo.bind('<<ComboboxSelected>>', lambda e: self.update_gel())
        
        # Botón actualizar
        ttk.Button(config_frame, text="Update Gel", 
                   command=self.update_gel).grid(row=4, column=0, columnspan=2, pady=10)
        
        # ===== GESTIÓN DE MUESTRAS =====
        samples_frame = ttk.LabelFrame(left_panel, text="Samples Management", padding=10)
        samples_frame.pack(fill='both', expand=True, padx=10, pady=5)
        
        # Entrada de nuevas muestras
        input_frame = ttk.Frame(samples_frame)
        input_frame.pack(fill='x', pady=5)
        
        ttk.Label(input_frame, text="Sample Name:").grid(row=0, column=0, sticky='w')
        self.sample_name = tk.StringVar()
        name_entry = ttk.Entry(input_frame, textvariable=self.sample_name, width=10)
        name_entry.grid(row=0, column=1, padx=5)
        
        ttk.Label(input_frame, text="Size (bp):").grid(row=0, column=2, sticky='w', padx=(10,0))
        self.sample_size = tk.StringVar()
        size_entry = ttk.Entry(input_frame, textvariable=self.sample_size, width=10)
        size_entry.grid(row=0, column=3, padx=5)
        
        ttk.Button(input_frame, text="Add Sample", 
                   command=self.add_sample).grid(row=0, column=4, padx=5)
        
        # Tabla de muestras
        columns = ('Sample Name', 'Size (bp)')
        self.tree = ttk.Treeview(samples_frame, columns=columns, show='headings', height=8)
        
        for col in columns:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=80)
        
        scrollbar = ttk.Scrollbar(samples_frame, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=scrollbar.set)
        
        self.tree.pack(side=tk.LEFT, fill='both', expand=True)
        scrollbar.pack(side=tk.RIGHT, fill='y')
        
        # Botones para manejar muestras
        button_frame = ttk.Frame(samples_frame)
        button_frame.pack(fill='x', pady=5)
        
        ttk.Button(button_frame, text="Delete Selected", 
                   command=self.remove_sample).pack(side=tk.LEFT, padx=2)
        ttk.Button(button_frame, text="Clear All", 
                   command=self.clear_samples).pack(side=tk.LEFT, padx=2)
        ttk.Button(button_frame, text="Load Example", 
                   command=self.load_example).pack(side=tk.LEFT, padx=2)
        
        # ===== PANEL DERECHO =====
        # Frame para el gráfico
        self.figure = Figure(figsize=(7, 6), dpi=100)
        self.ax = self.figure.add_subplot(111)
        
        self.canvas = FigureCanvasTkAgg(self.figure, right_panel)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Información adicional
        info_frame = ttk.Frame(right_panel)
        info_frame.pack(fill='x', pady=5)
        
        self.info_label = ttk.Label(info_frame, text="", font=('Arial', 9))
        self.info_label.pack()
    
    def add_sample(self):
        """Añade una muestra a la tabla"""
        name = self.sample_name.get().strip()
        size_str = self.sample_size.get().strip()
        
        if not name or not size_str:
            messagebox.showwarning("Input Error", "Please enter both name and size")
            return
        
        try:
            size = int(size_str)
            if size <= 0:
                raise ValueError("Size must be a positive integer")
        except ValueError:
            messagebox.showwarning("Input Error", "Size must be a positive integer")
            return
        
        # Verificar si el nombre ya existe
        for item in self.tree.get_children():
            if self.tree.item(item)['values'][0] == name:
                messagebox.showwarning("Duplicate Sample", f"A sample named '{name}' already exists")
                return
        
        # Añadir a la tabla
        self.tree.insert('', 'end', values=(name, size))
        self.samples.append({'name': name, 'size': size})
        
        # Limpiar campos
        self.sample_name.set('')
        self.sample_size.set('')
        
        # Actualizar gel
        self.update_gel()
    
    def remove_sample(self):
        """Elimina la muestra seleccionada"""
        selected = self.tree.selection()
        if not selected:
            messagebox.showwarning("Selection Error", "Please select a sample to delete")
            return
        
        # Eliminar de la tabla y de la lista
        for item in selected:
            values = self.tree.item(item)['values']
            self.tree.delete(item)
            self.samples = [s for s in self.samples if s['name'] != values[0]]
        
        self.update_gel()
    
    def clear_samples(self):
        """Limpia todas las muestras"""
        if messagebox.askyesno("Confirm Clear", "Delete all samples?"):
            for item in self.tree.get_children():
                self.tree.delete(item)
            self.samples = []
            self.update_gel()
    
    def load_example(self):
        """Carga un ejemplo de muestras"""
        self.clear_samples()
        
        example_samples = [
            ('Sample 1', 500),
            ('Sample 2', 1000),
            ('Sample 3', 750),
            ('Sample 4', 250),
            ('Sample 5', 1500),
            ('Sample 6', 300),
            ('Sample 7', 2000),
            ('Sample 8', 100)
        ]
        
        for name, size in example_samples:
            self.tree.insert('', 'end', values=(name, size))
            self.samples.append({'name': name, 'size': size})
        
        self.update_gel()
    
    def calculate_migration(self, fragment_size):
        """Calcula la distancia de migración para un fragmento"""
        # Factores de corrección por porcentaje de agarosa
        agarose_factors = {
            0.8: 1.0,
            1.0: 0.85,
            1.2: 0.75,
            1.5: 0.65,
            2.0: 0.55
        }
        
        agarose_factor = agarose_factors.get(self.gel_percentage.get(), 0.75)
        max_distance = 10.0
        voltage_factor = self.voltage.get() / 100
        time_factor = self.runtime.get() / 60
        
        # Fórmula de migración
        if fragment_size < 100:
            migration = max_distance * (1 - np.exp(-0.01 * fragment_size * agarose_factor))
        else:
            migration = max_distance * (1 - np.exp(-0.005 * fragment_size * agarose_factor))
        
        migration = migration * voltage_factor * time_factor
        return min(migration, max_distance * 1.5)
    
    def update_gel(self):
        """Actualiza la visualización del gel"""
        self.ax.clear()
        
        # Obtener datos
        ladder = self.ladders.get(self.selected_ladder.get(), self.ladders['1 kb'])
        samples = self.samples.copy()
        
        # Calcular migración para ladder
        ladder_data = []
        for size in ladder:
            migration = self.calculate_migration(size)
            if migration > 0:
                ladder_data.append({'size': size, 'migration': migration, 'intensity': 1.0})
        
        # Calcular migración para muestras
        samples_data = []
        for sample in samples:
            migration = self.calculate_migration(sample['size'])
            if migration > 0:
                samples_data.append({
                    'size': sample['size'],
                    'name': sample['name'],
                    'migration': migration,
                    'intensity': np.random.uniform(0.7, 1.0)
                })
        
        # Configurar el gel
        gel_width = 2.0
        max_migration = 15.0
        
        # Dibujar el gel
        def draw_lane(ax, x_pos, data, label, is_ladder=False, show_labels=True):
            # Fondo del carril
            ax.add_patch(Rectangle((x_pos - gel_width/2, 0), 
                                   gel_width, max_migration, 
                                   facecolor='lightgray', edgecolor='black', alpha=0.1))
            
            # Dibujar bandas
            for item in data:
                y_pos = max_migration - item['migration']
                if y_pos > 0.5 and y_pos < max_migration - 0.5:
                    # Ancho de la banda
                    band_width = gel_width * 0.8 * (0.5 + 0.5 * item['intensity'])
                    intensity = 0.3 + 0.7 * item['intensity']
                    color = (0, 0, 0, intensity)
                    
                    ax.add_patch(Rectangle((x_pos - band_width/2, y_pos - 0.08),
                                           band_width, 0.16,
                                           facecolor=color, edgecolor=color))
                    
                    # Etiquetas
                    if is_ladder and show_labels:
                        # Mostrar etiquetas para bandas principales de la ladder
                        if item['size'] in [10000, 5000, 3000, 2000, 1000, 500, 250, 100]:
                            ax.text(x_pos + gel_width/2 + 0.1, y_pos, 
                                   f"{item['size']}", fontsize=7, va='center')
                    elif not is_ladder and show_labels:
                        # Mostrar etiqueta con el tamaño
                        ax.text(x_pos + gel_width/2 + 0.1, y_pos, 
                               f"{item['size']}", fontsize=7, va='center')
            
            # Etiqueta del carril
            ax.text(x_pos, max_migration + 0.5, label, ha='center', fontsize=9)
        
        # Dibujar carril de ladder
        draw_lane(self.ax, 0, ladder_data, f"Ladder\n{self.selected_ladder.get()}", 
                  is_ladder=True, show_labels=True)
        
        # Dibujar carriles para muestras
        if samples_data:
            # Ordenar por tamaño para visualización
            sorted_samples = sorted(samples_data, key=lambda x: x['size'])
            
            for i, sample in enumerate(sorted_samples):
                x_pos = (i + 1) * 2.5
                draw_lane(self.ax, x_pos, [sample], sample['name'], 
                          is_ladder=False, show_labels=True)
        
        # Configurar el gráfico
        max_lanes = max(1, len(samples_data))
        self.ax.set_xlim(-2, max(10, (max_lanes + 1) * 2.5))
        self.ax.set_ylim(-1, max_migration + 2)
        self.ax.set_xlabel('Lanes')
        self.ax.set_ylabel('Migration Distance (cm)')
        self.ax.set_title(f'Agarose Gel {self.gel_percentage.get()}% - {self.voltage.get()}V, {self.runtime.get()}min')
        self.ax.grid(True, alpha=0.2)
        self.ax.invert_yaxis()
        
        # Actualizar información
        info_text = f"Samples: {len(samples_data)} | Ladder: {self.selected_ladder.get()} | "
        if samples_data:
            sizes = [str(s['size']) for s in samples_data[:5]]
            info_text += f"Sizes: {', '.join(sizes)}"
            if len(samples_data) > 5:
                info_text += ", ..."
        else:
            info_text += "No samples loaded"
        self.info_label.config(text=info_text)
        
        # Redibujar
        self.canvas.draw()
    
    def export_gel(self):
        """Exporta el gel como imagen"""
        from tkinter import filedialog
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("All files", "*.*")]
        )
        if filename:
            self.figure.savefig(filename, dpi=300, bbox_inches='tight')
            messagebox.showinfo("Export Successful", f"Image saved to: {filename}")

# Función principal
def main():
    root = tk.Tk()
    app = PCRGelApp(root)
    
    # Agregar menú de exportación
    menubar = tk.Menu(root)
    file_menu = tk.Menu(menubar, tearoff=0)
    file_menu.add_command(label="Export Gel Image", command=app.export_gel)
    file_menu.add_separator()
    file_menu.add_command(label="Exit", command=root.quit)
    menubar.add_cascade(label="File", menu=file_menu)
    root.config(menu=menubar)
    
    root.mainloop()

if __name__ == "__main__":
    main()