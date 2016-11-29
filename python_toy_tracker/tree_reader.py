import ROOT
import getopt
import sys

input_file = "/home/john/smear_missed_straws_5_track_10_fit.root"

# Get arguments and options from console
try:
    opts, args = getopt.getopt(argv, "hi", ["input_file="])
except getopt.GetoptError:
    print "tree_reader.py -i <input_file>"
    sys.exit(2)

# Set parameters according to arguments
for opt, arg in opts:
    if opt = "-h":
        print "tree_reader.py -i <input_file>"
        sys.exit()
    elif opt in ("-i", "--input_file"):
        input_file = arg


# Get input file
f = ROOT.TFile.Open(input_file, "read")

# Get tree of straw and scintillator strikes
tree = f.Get("event_tree")

tree.Print()
print ""
print "Full Tree Entries:", tree.GetEntries()

# Get entries from first fit
cut_tree = tree.CopyTree('fitNum==0')

print "Cut Tree Entries:", cut_tree.GetEntries()

# Loop across entries for first fit, and print
for i in xrange(cut_tree.GetEntries()):

    cut_tree.GetEntry(i)
 
    print "Cut_Tree Entry:", str(i)
    print "Indexes:", str(cut_tree.fitNum), str(cut_tree.eventNum)
    print "Alignments:,", str(cut_tree.trueAlignment), str(cut_tree.fittedAlignment) 
    print "True Track Pars:", cut_tree.trueTrackGrad, cut_tree.trueTrackInt
    print "Fitted Track Pars:", cut_tree.fittedTrackGrad, cut_tree.fittedTrackInt
    print "True Hit Dist:", cut_tree.trueHitDistance
    print "Fitted Hit Dist:", cut_tree.fittedHitDistance
    print "Module Num:", cut_tree.moduleNum
    print "Plane Num:", cut_tree.planeNum
    print "Layer Num:", cut_tree.layerNum
    print "Wire Num:", cut_tree.wireNum
    print ""
