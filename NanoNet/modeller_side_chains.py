from modeller import *
from modeller.automodel import *
import sys

def make_alignment_file(pdb_name, sequence):

    ali_file = open("pdb_seq.ali", "w") 
    ali_file.write(">P1;pdb_seq\n")
    ali_file.write("sequence:pdb_seq:::::::0.00: 0.00\n")
    ali_file.write(sequence + "*")
    ali_file.close()
    
    pdb_file = pdb_name
    env = environ()
    aln = alignment(env)
    mdl = model(env, file=pdb_file)
    aln.append_model(mdl, align_codes=pdb_file, atom_files=pdb_file)
    aln.append(file=f"pdb_seq.ali", align_codes="pdb_seq")
    aln.align2d()
    aln.write(file="alignment_for_full_atoms.ali", alignment_format='PIR')


class MyLoopModel(automodel):
    def special_patches(self, aln):
        # Rename heavy chain and renumber the residues
        self.rename_segments(segment_ids=['H'],renumber_residues=[1])

def relax_pdb(pdb_name, sequence):

    make_alignment_file(pdb_name, sequence)

    pdb_file = pdb_name
    log.verbose()
    env = environ()

    # directories for input atom files
    env.io.atom_files_directory = ['.', '../atom_files']

    a = MyLoopModel(env, alnfile='alignment_for_full_atoms.ali', knowns=pdb_file, sequence="pdb_seq")
    a.starting_model = 1
    a.ending_model = 1
    a.make()




relax_pdb(pdb_name=sys.argv[1], sequence=sys.argv[2])
