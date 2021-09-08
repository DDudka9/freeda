"""

This module will have to be incorporated into __main__ (freeda_pipeline.py) cose it calls it

"""

from tkinter import *
from tkinter import ttk
import os
import re

wdir = os.getcwd() + "/"

error_message = StringVar()
message = "Invalid gene name (follow pattern: Cenpo)"

def make_proteins_file(*args):
    with open(wdir + "proteins_new.txt", "a") as f:
        f.write("\n" + value)


def check_gene_name(gene_name, op):
    """Checks if user provided a valid gene name"""
    error_message.set("")
    valid = re.match(r"^[A-Z]{1}([A-Za-z0-9]+$)", gene_name) is not None
    btn.state(["!disabled"] if valid else ["disabled"])
    if op=="key":
        ok_so_far = re.match


root = Tk()
root.title("FREEDA - Finder of Rapidly Evolving Exons in De novo Assemblies")
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

#s = ttk.Style()
#s.configure("Danger.TFrame", background="red", borderwidth=5, relief="raised")

#mainframe = ttk.Frame(root, width=200, height=299, style="Danger.TFrame", padding="5 5 5 5")
mainframe = ttk.Frame(root, width=200, height=299, padding="5 5 5 5")
mainframe.grid(column=5, row=20, sticky=(N, W, E, S))
#mainframe["borderwidth"] = 2
#mainframe["relief"] = "sunken"

# this will be a gene name provided by the user
#protein = StringVar()
#protein_to_analyze = ttk.Label(mainframe, width=20, text="Gene name :", justify="left", textvariable=protein)
#protein_to_analyze.grid(column=1, row=1, sticky=(W, E))
# set the new protein name following user input
#protein.set(user_gene_name)


# user chooses which blast threshold should be used

t = IntVar()
ttk.Label(mainframe, text="Blast threshold").grid(column=1, row=1, sticky=(E))
low = ttk.Radiobutton(mainframe, text="30%", variable=t, value="low").grid(column=2, row=2, sticky=(W))
medium = ttk.Radiobutton(mainframe, text="50%", variable=t, value="medium").grid(column=2, row=3, sticky=(W))
high = ttk.Radiobutton(mainframe, text="70%", variable=t, value="high").grid(column=2, row=4, sticky=(W))

# user chooses which clade to analyze

clade = StringVar()
ttk.Label(mainframe, text="Clade").grid(column=1, row=5, sticky=(E))
rodents = ttk.Radiobutton(mainframe, text="rodents", variable=clade, value="rodents").grid(column=2, row=6, sticky=(W))
primates = ttk.Radiobutton(mainframe, text="primates", variable=clade, value="primates").grid(column=2, row=7, sticky=(W))

# user inputs up to 5 gene names to analyse

gene_name1 = StringVar()
ttk.Label(mainframe, text="Gene name").grid(column=1, row=8, sticky=(E))
protein1 = ttk.Entry(mainframe, textvariable=gene_name1).grid(column=2, row=8, sticky=(W))

gene_name2 = StringVar()
ttk.Label(mainframe, text="Gene name").grid(column=1, row=9, sticky=(E))
protein2 = ttk.Entry(mainframe, textvariable=gene_name2).grid(column=2, row=9, sticky=(W))

gene_name3 = StringVar()
ttk.Label(mainframe, text="Gene name").grid(column=1, row=10, sticky=(E))
protein3 = ttk.Entry(mainframe, textvariable=gene_name3).grid(column=2, row=10, sticky=(W))

gene_name4 = StringVar()
ttk.Label(mainframe, text="Gene name").grid(column=1, row=11, sticky=(E))
protein4 = ttk.Entry(mainframe, textvariable=gene_name4).grid(column=2, row=11, sticky=(W))

gene_name5 = StringVar()
ttk.Label(mainframe, text="Gene name").grid(column=1, row=12, sticky=(E))
protein5 = ttk.Entry(mainframe, textvariable=gene_name5).grid(column=2, row=12, sticky=(W))

button = ttk.Button(mainframe, text="Analyze", default="active").grid(column=5, row=20, sticky=(W)) #, command=freeda.freeda_pipeline)

#threshold = StringVar()
#threshold_label = ttk.Label(mainframe, textvariable=threshold).grid(column=1, row=1, sticky=(E))
#threshold_entry = ttk.Entry(mainframe, textvariable=threshold, width=2).grid(column=1, row=1, sticky=(W))


# make button to run the pipeline
#button = ttk.Button(mainframe, text="Analyze", default="active", command=freeda.freeda_pipeline)
# will execute script attached to the button when left mouse clicked
#root.bind("<ButtonPress-1>", lambda e: button.invoke())

root.mainloop()

"""

root = Tk()
l = ttk.Label(root, text="Starting...")
l.grid()
l.bind("<Enter>", lambda e: l.configure(text="Moved mouse inside"))
l.bind("<Leave>", lambda e: l.configure(text="Moved mouse outside"))
l.bind("<ButtonPress-1>", lambda e: l.configure(text="Clicked left mouse button"))
l.bind("<ButtonPress-2>", lambda e: l.configure(text="Clicked right mouse button"))
l.bind("<Double-1>", lambda e: l.configure(text="Double clicked"))
l.bind("<Motion>", lambda e: l.configure(text="Mouse moved to %d, %d" % (e.x, e.y)))

root.mainloop()


def calculate(*args):
    try:
        value = float(feet.get())
        meters.set(int(0.3048 * value * 10000.0 + 0.5)/10000.0)
    except ValueError:
        pass

root = Tk()
root.title("Feet to Meters")

mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

feet = StringVar()
feet_entry = ttk.Entry(mainframe, width=7, textvariable=feet)
feet_entry.grid(column=2, row=1, sticky=(W, E))

meters = StringVar()
ttk.Label(mainframe, textvariable=meters).grid(column=2, row=2, sticky=(W, E))

ttk.Button(mainframe, text="Calculate", command=calculate).grid(column=3, row=3, sticky=(W))

ttk.Label(mainframe, text="feet").grid(column=3, row=1, sticky=(W))
ttk.Label(mainframe, text="is equivalent to").grid(column=1, row=2, sticky=(E))
ttk.Label(mainframe, text="meters").grid(column=3, row=2, sticky=(W))

for child in mainframe.winfo_children():
    child.grid_configure(padx=5, pady=5)

feet_entry.focus()
root.bind("<Return>", calculate)

root.mainloop()

"""