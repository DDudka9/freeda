"""

Modified from Roger Edgar (robert@drive5.com) work-around to -stable bug option for MUSCLE aligner
Keeps aligned sequences in order
Taken from: https://www.drive5.com/muscle/manual/stable.html

"""

#import sys

def stable_order(input, output):

	#InputFileName = sys.argv[1]
	#AlnFileName = sys.argv[2]

	InLabels, InSeqs = ReadSeqs2(input)
	AlnSeqs = ReadSeqs(output)

	for Label in InLabels:
		if Label not in AlnSeqs.keys():
			print("Not found in alignment: " + Label)
		print(">" + Label)
		print(AlnSeqs[Label])

#def Die(s):
	#print(sys.stderr)
	#print(sys.stderr, sys.argv)
	#print(sys.stderr, "**ERROR**", s)
	#sys.exit(1)


def ReadSeqs(output):
	Seqs = {}
	Id = ""
	File = open(output)
	while 1:
		Line = File.readline()
		if len(Line) == 0:
			return Seqs
		if len(Line) == 0:
			continue
		if Line[0] == ">":
			Id = Line[1:][:-1]
			Seqs[Id] = ""
		else:
			if Id == "":
				print("FASTA file '%s' does not start with '>'" % output)
			Seqs[Id] = Seqs[Id] + Line[:-1]

	# write into file
	with open(output, "w") as f:
		for id, seq in Seqs:
			f.write(id)
			f.write(seq)


def ReadSeqs2(input):
	Seqs = []
	Labels = []
	File = open(input)
	while 1:
		Line = File.readline()
		if len(Line) == 0:
			return Labels, Seqs
		Line = Line.strip()
		if len(Line) == 0:
			continue
		if Line[0] == ">":
			Id = Line[1:]
			Labels.append(Id)
			Seqs.append("")
		else:
			i = len(Seqs)-1
			Seqs[i] = Seqs[i] + Line

