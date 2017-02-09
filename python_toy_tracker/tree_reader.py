import ROOT
import getopt
import sys

input_file = "/home/john/smear_missed_straws_5_track_10_fit.root"

# Get arguments and options from console


argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "h:i:", ["help", "input_file="])
except getopt.GetoptError:
    print "tree_reader.py -i <input_file>"
    sys.exit(2)

# Set parameters according to arguments
for opt, arg in opts:

    if opt in ("-h", "--help"):
        print "tree_reader.py -i <input_file>"
        sys.exit()
    elif opt in ("-i", "--input_file"):
        input_file = arg


print input_file

# Get input file
f = ROOT.TFile.Open(input_file, "read")

# Get tree of straw hit events
tree = f.Get("event_tree")


print ""
print "Full Tree Entries:", tree.GetEntries()

# Get entries from first fit
#cut_tree = tree.CopyTree('fitNum==0')

#print "Cut Tree Entries:", cut_tree.GetEntries()

# Loop across entries for first fit, and print
for i in xrange(tree.GetEntries()):

    tree.GetEntry(i)
 
    print "Tree Entry:", str(i)
    print "Indexes:", str(tree.fitNum), str(tree.eventNum)
    print "Alignments:,", str(tree.trueAlignment), str(tree.fittedAlignment) 
    print "True Track Pars:", tree.trueTrackGrad, tree.trueTrackInt
    print "Fitted Track Pars:", tree.fittedTrackGrad, tree.fittedTrackInt
    print "True Hit Dist:", tree.trueHitDistance
    print "Fitted Hit Dist:", tree.fittedHitDistance
    print "Module Num:", tree.moduleNum
    print "Plane Num:", tree.planeNum
    print "Layer Num:", tree.layerNum
    print "Wire Num:", tree.wireNum
    print ""
