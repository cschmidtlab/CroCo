import sys, os

croco_script = sys.argv[0]
# check if programme was called via symlink
if os.path.islink(croco_script):
    croco_script = os.readlink(croco_script)
# dir is the directory above the bin-dir
croco_dir = os.path.abspath(os.path.join(os.path.dirname(croco_script), '..'))

sys.path.append(os.path.abspath(os.path.join('..', croco_dir)))
    
import tkinter as tk
from tkinter import ttk

import croco
from pandas import read_csv

class CroCoMainWindow():
    def __init__(self):
        self.win = tk.Tk()
        self.win.title("The CroCo Crosslink converter")      
        self.availReads = {'pLink1': croco.pLink1.Read,
                           'pLink2': croco.pLink2.Read,
                           'Kojak': croco.Kojak.Read,
                           'xQuest': croco.xQuest.Read,
                           'xTable': read_csv}
        self.availWrites =  {'xTable': croco.xTable.Write,
                             'xVis': croco.xVis.Write,
                             'xiNet': croco.xiNET.Write,
                             'DynamXL': croco.DynamXL.Write,
                             'xWalk': croco.xWalk.Write}
        self.createWidgets()

    def createWidgets(self):
        ## Input Frame
        self.input_frame = tk.LabelFrame(self.win, text='Input')
        self.input_frame.grid(column=0, row=0)
        
        # Select input format
        ttk.Label(self.input_frame, text='Input file format').grid(column=0, row=0)
        self.input_format = tk.StringVar()
        self.input_dropdown = ttk.Combobox(self.input_frame, textvariable=self.input_format, state='readonly')
        self.input_dropdown['values'] = list(self.availReads.keys())
        self.input_dropdown.grid(column=1, row=0)
        
        # Select output format
        ttk.Label(self.input_frame, text='Output file format').grid(column=2, row=0)
        self.output_format = tk.StringVar()
        self.input_dropdown = ttk.Combobox(self.input_frame, textvariable=self.output_format, state='readonly')
        self.input_dropdown['values'] = list(self.availWrites.keys())
        self.input_dropdown.grid(column=3, row=0)
        
        self.load_btn = ttk.Button(self.input_frame, 'Load file(s)')
        
        ## Control Frame
        control_frame = tk.Frame(self.win)
        control_frame.grid(column=0, row=1)

main = CroCoMainWindow()

# Start GUI
main.win.mainloop()