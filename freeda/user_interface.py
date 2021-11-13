"""

This module will have to be incorporated into __main__ (freeda_pipeline.py) cose it calls it

"""

from tkinter import *
from tkinter import ttk
from freeda import input_extractor
import os
import re
import pyensembl


# Define which operating system is used -> os.uname().sysname
# Use spinbox for gene names? -> it will contain all possible pyensembl genes for the release
# Use scrollbar?


def check_gene_name(gene_name, op):
    """Checks if user provided a valid gene name"""
    error_message1.set("")
    # accept only entry starting with one capital letter followed by small or big letters or numbers
    valid = re.match(r"^[A-Z]{1}([A-Za-z0-9]+$)", gene_name) is not None
    # button can be clicked only if gene names are valid
    button.state(["!disabled"] if valid else ["disabled"])
    # keystroke validation
    if op == "key":
        ok_so_far = re.match(r"[A-Za-z0-9]+$", gene_name) is not None
        if not ok_so_far:
            error_message1.set(message1)
        return ok_so_far
    elif op == "focusout":
        if not valid:
            error_message1.set(message1)
    return valid


def check_functional_residues(residue, op):
    """Checks if user provided a valid residue number"""
    error_message2.set("")
    # accept only entry starting with one capital letter followed by small or big letters or numbers
    valid = re.match(r"^[1-9]{1}([0-9]{0,3}$)", residue) is not None  # max 4 digits, first one not a zero
    # button can be clicked only if gene names are valid
    button.state(["!disabled"] if valid else ["disabled"])
    # keystroke validation
    if op == "key":
        ok_so_far = re.match(r"[0-9]+$", residue) is not None
        if not ok_so_far:
            error_message2.set(message2)
        return ok_so_far
    elif op == "focusout":
        if not valid:
            error_message2.set(message2)
    return valid


wdir = os.getcwd() + "/"


# this will be a gene name provided by the user
#protein = StringVar()
#protein_to_analyze = ttk.Label(mainframe, width=20, text="Gene name :", justify="left", textvariable=protein)
#protein_to_analyze.grid(column=1, row=1, sticky=(W, E))
# set the new protein name following user input
#protein.set(user_gene_name)


#ttk.Label(mainframe, text="Genes to analyze").grid()


# set up the main window
root = Tk()
root.title("FREEDA - Finder of Rapidly Evolving Exons in De novo Assemblies")
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

error_message1 = StringVar()
message1 = "Invalid gene name (follow pattern: Cenpo for rodents and : CENPO for primates)"
error_message2 = StringVar()
message2 = "Invalid residue number (follow pattern: 100)"

# create a mainframe inside the parent frame
mainframe = ttk.Frame(root, width=200, height=400, padding="5 5 5 5")
mainframe.grid(sticky=(N, W, E, S))

# configure first two columns of the mainframe to behave uniformly
mainframe.columnconfigure(0, weight=1, uniform="group1")
mainframe.columnconfigure(1, weight=1, uniform="group1")
mainframe.rowconfigure(0, weight=1)

# create input frame
input_frame = ttk.Frame(mainframe, relief="ridge", padding="5 5 5 5")
input_frame.grid(column=0, row=0, sticky=(N, W, E, S), padx=5, pady=5)
# let all columns resize
input_frame.columnconfigure(0, weight=3, uniform="group1")
input_frame.columnconfigure((1, 2), weight=1, uniform="group1")
input_frame.columnconfigure(3, weight=2, uniform="group1")


# create output frame
output_frame = ttk.Frame(mainframe, relief="ridge", padding="5 5 5 5")
output_frame.grid(column=1, row=0, sticky=(N, W, E, S), padx=5, pady=5)
# let all columns resize
output_frame.columnconfigure(0, weight=3, uniform="group1")
output_frame.columnconfigure((1, 2), weight=1, uniform="group1")
output_frame.columnconfigure(3, weight=2, uniform="group1")

# create and grid the logging window
logging_window = Text(output_frame, state="disabled", wrap="none")
logging_window.grid(column=0, row=1, columnspan=4, sticky=(N, W, E, S))
ttk.Label(output_frame, text="Events window (logged to 'FREEDA*.log'").grid(column=0, row=0, columnspan=4, sticky=(W))

#mainframe["borderwidth"] = 2
#mainframe["relief"] = "sunken"

# user chooses which clade to analyze

clade = StringVar()
ttk.Label(input_frame, text="Clade").grid(column=0, row=0, pady=5, sticky=(W))
rodents = ttk.Radiobutton(input_frame, text="Rodents", variable=clade, value="Rodents")
rodents.grid(column=0, row=1, sticky=(W))
primates = ttk.Radiobutton(input_frame, text="Primates", variable=clade, value="Primates")
primates.grid(column=0, row=2, sticky=(W))

# user chooses which blast threshold should be used

threshold = IntVar()
ttk.Label(input_frame, text="Blast search").grid(column=0, row=4, pady=5, sticky=(W))
deep = ttk.Radiobutton(input_frame, text="Deep", variable=threshold,
                       value="deep").grid(column=0, row=5, sticky=(W))
medium = ttk.Radiobutton(input_frame, text="Medium (recommended)", variable=threshold,
                         value="medium").grid(column=0, row=6, sticky=(W))
shallow = ttk.Radiobutton(input_frame, text="Shallow", variable=threshold,
                          value="shallow").grid(column=0, row=7, sticky=(W))


check_gene_name_wrapper = (input_frame.register(check_gene_name), "%P", "%V")
check_functional_residues = (input_frame.register(check_functional_residues), "%P", "%V")

# user inputs up to 3 gene names to analyse

release = 104
species = "mus musculus"
ensembl = pyensembl.EnsemblRelease(release, species)
ensembl.download()  # this is suppose to bypass installing the release from outside python
ensembl.index()  # this is suppose to bypass installing the release from outside python



#update(all_genes)

ttk.Label(input_frame, text="Indicate functional residues").grid(column=1, row=8, columnspan=3, sticky=(N))
ttk.Label(input_frame, text="start").grid(column=1, row=9)
ttk.Label(input_frame, text="end").grid(column=2, row=9)
ttk.Label(input_frame, text="label").grid(column=3, row=9)

# GENE 1

gene1_frame = ttk.Frame(input_frame, relief="ridge", padding="5 5 5 5")
gene1_frame.grid(column=0, row=11, columnspan=4, sticky=(N, W, E, S), padx=5, pady=5)
# let all columns resize
gene1_frame.columnconfigure(0, weight=3, minsize=50, uniform="group1")
gene1_frame.columnconfigure((1, 2), weight=1, minsize=50, uniform="group1")
gene1_frame.columnconfigure(3, weight=2, minsize=50, uniform="group1")


gene_name1 = StringVar()
ttk.Label(gene1_frame, text="Gene name").grid(column=0, row=11, padx=6, sticky=(W))
#gene1_list = Tk.ListBox(mainframe, height=10, width=10)
#gene1_list.pack(pady=40)
#all_genes = input_extractor.make_gene_list(ensembl)
protein1 = ttk.Entry(gene1_frame, textvariable=gene_name1, validate="all", validatecommand=check_gene_name_wrapper)
protein1.grid(column=0, row=12, padx=5, pady=2, sticky=(W))
dup1 = StringVar()
ttk.Checkbutton(gene1_frame, text="Duplication expected", #command=duplication_expected,
                              variable=dup1, onvalue="on", offvalue="off").grid(column=0, row=13, padx=6, sticky=(W))
s11_start = StringVar()
s11_end = StringVar()
s11_label = StringVar()
site11_start = ttk.Entry(gene1_frame, textvariable=s11_start, validate="all", validatecommand=check_functional_residues)
site11_start.grid(column=1, row=11, padx=5, pady=2, sticky=(W))
site11_end = ttk.Entry(gene1_frame, textvariable=s11_end, validate="all", validatecommand=check_functional_residues)
site11_end.grid(column=2, row=11, padx=5, pady=2, sticky=(W))
site11_label = ttk.Entry(gene1_frame, textvariable=s11_label, validate="all")
site11_label.grid(column=3, row=11, padx=5, pady=2, sticky=(W))

s12_start = StringVar()
s12_end = StringVar()
s12_label = StringVar()
site12_start = ttk.Entry(gene1_frame, textvariable=s12_start, validate="all", validatecommand=check_functional_residues)
site12_start.grid(column=1, row=12, padx=5, pady=2, sticky=(W))
site12_end = ttk.Entry(gene1_frame, textvariable=s12_end, validate="all", validatecommand=check_functional_residues)
site12_end.grid(column=2, row=12, padx=5, pady=2, sticky=(W))
site12_label = ttk.Entry(gene1_frame, textvariable=s12_label, validate="all")
site12_label.grid(column=3, row=12, padx=5, pady=2, sticky=(W))

s13_start = StringVar()
s13_end = StringVar()
s13_label = StringVar()
site13_start = ttk.Entry(gene1_frame, textvariable=s13_start, validate="all", validatecommand=check_functional_residues)
site13_start.grid(column=1, row=13, padx=5, pady=2, sticky=(W))
site13_end = ttk.Entry(gene1_frame, textvariable=s13_end, validate="all", validatecommand=check_functional_residues)
site13_end.grid(column=2, row=13, padx=5, pady=2, sticky=(W))
site13_label = ttk.Entry(gene1_frame, textvariable=s13_label, validate="all")
site13_label.grid(column=3, row=13, padx=5, pady=2, sticky=(W))

# GENE 2

gene2_frame = ttk.Frame(input_frame, relief="ridge", padding="5 5 5 5")
gene2_frame.grid(column=0, row=14, columnspan=4, sticky=(N, W, E, S), padx=5, pady=5)
# let all columns resize
gene2_frame.columnconfigure(0, weight=3, uniform="group1")
gene2_frame.columnconfigure((1, 2), weight=1, uniform="group1")
gene2_frame.columnconfigure(3, weight=2, uniform="group1")

gene_name2 = StringVar()
ttk.Label(gene2_frame, text="Gene name").grid(column=0, row=14, padx=6, sticky=(W))
protein2 = ttk.Entry(gene2_frame, textvariable=gene_name2, validate="all", validatecommand=check_gene_name_wrapper)
protein2.grid(column=0, row=15, padx=5, pady=5, sticky=(W))
dup2 = StringVar()
ttk.Checkbutton(gene2_frame, text="Duplication expected", #command=duplication_expected,
                              variable=dup2, onvalue="on", offvalue="off").grid(column=0, row=16, padx=6, sticky=(W))

s21_start = StringVar()
s21_end = StringVar()
s21_label = StringVar()
site21_start = ttk.Entry(gene2_frame, textvariable=s21_start, validate="all", validatecommand=check_functional_residues)
site21_start.grid(column=1, row=14, padx=5, pady=2, sticky=(W))
site21_end = ttk.Entry(gene2_frame, textvariable=s21_end, validate="all", validatecommand=check_functional_residues)
site21_end.grid(column=2, row=14, padx=5, pady=2, sticky=(W))
site21_label = ttk.Entry(gene2_frame, textvariable=s21_label, validate="all")
site21_label.grid(column=3, row=14, padx=5, pady=2, sticky=(W))

s22_start = StringVar()
s22_end = StringVar()
s22_label = StringVar()
site22_start = ttk.Entry(gene2_frame, textvariable=s22_start, validate="all", validatecommand=check_functional_residues)
site22_start.grid(column=1, row=15, padx=5, pady=2, sticky=(W))
site22_end = ttk.Entry(gene2_frame, textvariable=s22_end, validate="all", validatecommand=check_functional_residues)
site22_end.grid(column=2, row=15, padx=5, pady=2, sticky=(W))
site22_label = ttk.Entry(gene2_frame, textvariable=s22_label, validate="all")
site22_label.grid(column=3, row=15, padx=5, pady=2, sticky=(W))

s23_start = StringVar()
s23_end = StringVar()
s23_label = StringVar()
site23_start = ttk.Entry(gene2_frame, textvariable=s23_start, validate="all", validatecommand=check_functional_residues)
site23_start.grid(column=1, row=16, padx=5, pady=2, sticky=(W))
site23_end = ttk.Entry(gene2_frame, textvariable=s23_end, validate="all", validatecommand=check_functional_residues)
site23_end.grid(column=2, row=16, padx=5, pady=2, sticky=(W))
site23_label = ttk.Entry(gene2_frame, textvariable=s23_label, validate="all")
site23_label.grid(column=3, row=16, padx=5, pady=2, sticky=(W))

# GENE 3

gene3_frame = ttk.Frame(input_frame, relief="ridge", padding="5 5 5 5")
gene3_frame.grid(column=0, row=17, columnspan=4, sticky=(N, W, E, S), padx=5, pady=5)
# let all columns resize
gene3_frame.columnconfigure(0, weight=3, uniform="group1")
gene3_frame.columnconfigure((1, 2), weight=1, uniform="group1")
gene3_frame.columnconfigure(3, weight=2, uniform="group1")

gene_name3 = StringVar()
ttk.Label(gene3_frame, text="Gene name").grid(column=0, row=17, padx=6, sticky=(W))
protein3 = ttk.Entry(gene3_frame, textvariable=gene_name3, validate="all", validatecommand=check_gene_name_wrapper)
protein3.grid(column=0, row=18, padx=5, pady=5, sticky=(W))
dup3 = StringVar()
ttk.Checkbutton(gene3_frame, text="Duplication expected", #command=duplication_expected,
                              variable=dup3, onvalue="on", offvalue="off").grid(column=0, row=19, padx=6, sticky=(W))

s31_start = StringVar()
s31_end = StringVar()
s31_label = StringVar()
site31_start = ttk.Entry(gene3_frame, textvariable=s31_start, validate="all", validatecommand=check_functional_residues)
site31_start.grid(column=1, row=17, padx=5, pady=2, sticky=(W))
site31_end = ttk.Entry(gene3_frame, textvariable=s31_end, validate="all", validatecommand=check_functional_residues)
site31_end.grid(column=2, row=17, padx=5, pady=2, sticky=(W))
site31_label = ttk.Entry(gene3_frame, textvariable=s31_label, validate="all")
site31_label.grid(column=3, row=17, padx=5, pady=2, sticky=(W))

s32_start = StringVar()
s32_end = StringVar()
s32_label = StringVar()
site32_start = ttk.Entry(gene3_frame, textvariable=s32_start, validate="all", validatecommand=check_functional_residues)
site32_start.grid(column=1, row=18, padx=5, pady=2, sticky=(W))
site32_end = ttk.Entry(gene3_frame, textvariable=s32_end, validate="all", validatecommand=check_functional_residues)
site32_end.grid(column=2, row=18, padx=5, pady=2, sticky=(W))
site32_label = ttk.Entry(gene3_frame, textvariable=s32_label, validate="all")
site32_label.grid(column=3, row=18, padx=5, pady=2, sticky=(W))

s33_start = StringVar()
s33_end = StringVar()
s33_label = StringVar()
site33_start = ttk.Entry(gene3_frame, textvariable=s33_start, validate="all", validatecommand=check_functional_residues)
site33_start.grid(column=1, row=19, padx=5, pady=2, sticky=(W))
site33_end = ttk.Entry(gene3_frame, textvariable=s33_end, validate="all", validatecommand=check_functional_residues)
site33_end.grid(column=2, row=19, padx=5, pady=2, sticky=(W))
site33_label = ttk.Entry(gene3_frame, textvariable=s33_label, validate="all")
site33_label.grid(column=3, row=19, padx=5, pady=2, sticky=(W))

# search for working directory (Data folder)

wdir_label = StringVar()
site11_start = ttk.Entry(gene1_frame, textvariable=s11_start, validate="all", validatecommand=check_functional_residues)
site11_start.grid(column=1, row=11, padx=5, pady=2, sticky=(W))



button = ttk.Button(mainframe, text="Analyze") # default="active"
button.grid(column=5, row=20, padx=5, pady=5, sticky=(W)) #, command=freeda.freeda_pipeline)
button.state(["disabled"])

error_label1 = ttk.Label(mainframe, font="TkSmallCaptionFont", foreground="red", textvariable=error_message1)
error_label1.grid(column=2, row=8, padx=5, pady=5, sticky="w")

error_label2 = ttk.Label(mainframe, font="TkSmallCaptionFont", foreground="red", textvariable=error_message2)
error_label2.grid(column=2, row=9, padx=5, pady=5, sticky="w")

#ttk.Label(mainframe, text="Gene name").grid(column=0, row=11, sticky=(W))

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