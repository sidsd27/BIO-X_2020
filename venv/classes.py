
class guide_RNA:
    def __init__(self, crRNA, PAM, d_insert, hUS, hDS, strand):
        self.PAM = PAM
        self.d_insert = d_insert
        self.crRNA = crRNA
        self.hUS = hUS
        self.hDS = hDS
        self.strand = strand


    def print_features(self, mutated):
        print("Sequence", self.crRNA, "PAM", self.PAM, "Cut site", self.d_insert)
        print("57 base pair Homology Upstream:", self.hUS, "\n")
        if mutated == 1:
         print("57 base pair Homology Downstream (mutated):" + str(self.hDS[0]) + "GT" + str(self.hDS[3:]))
         print("57 base pair Homology Downstream:",self.hDS,"\n")

    def print(self,x):
        print(self.x)

class guide_data:
    def __init__(self, sequence, target_found, distance, activity, off_target, common_exon, hUS, hDS):
        self.sequence = sequence
        self.distance = distance
        self.activity = activity
        self.off_target = off_target
        self.common_exon = common_exon
        self.target_found = target_found
        self.hUS = hUS
        self.hDS = hDS


    def print_guide_data(self):
        print("Sequence:", self.sequence)
        print("Target found:", self.target_found)
        print("Distance from insertion site:", self.distance)
        print("Low off-target activity:", self.off_target)
        print("High activity:", self.activity)
        print("LHA:", self.hUS)
        print("RHA:", self.hDS)
        print("Targets common exon in all transcripts:", self.common_exon, "\n")
