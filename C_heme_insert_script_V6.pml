###PyMOL c-type cytochrome heme insertion script for AlphaFold2 models written by Benjamin W. Nash at Univeristy of East Anglia School of Biological Sciences, Norwich, United Kingdom, NR4 7TJ, in 2022. This work was made possible through BBSRC grant BB/T008717/1.
###Should you have problems using the script and wish to contact me for help, please get in touch with me at Benjamin.Nash5@outlook.com.

########## if you wish for the hemes to be labelled under a specific chain name, please change the word Heme_Chain in the following line, to your chain word of choice ######
G = "Heme_Chain"

sele ch_motifs, pepseq CH, 
#make selection called ch_motifs by selecting all CH motifs

sele cxxch_motifs, ch_motifs or (r. CYS and (ch_motifs xt. 20))
#expand selection to include the first cysteine

python
a = cmd.count_atoms("cxxch_motifs and name ca")
number_of_hemes = int(a / 3)
python end
print number_of_hemes
#counts calpha carbons in cxxch motifs and calculates number of hemes (3 calpha per CXXCH)

sele cxxch_histidines, cxxch_motifs (and r. HIS (and n. CA))
#creates selection of cxxch histidine calpha carbons for iterating the alignment

fetch 6HR0
#get pdb STC

time.sleep(5)

select remove_from_heme_to_add, /6HR0 and (not resi 101 (or resi 58) (or resi 61) (or resi 62))
#make selection of STC to remove everything that isn't heme 1

alter /6HR0/B/A/HEC`101,segi=''
#set the STC heme to no segment, will be helpful for later
alter /6HR0//A/HEC`101,resi='3000'
#set the STC heme to residue 3000, will be helpful for later

remove (remove_from_heme_to_add)
#actually remove everything that isn't heme + cxxch

delete remove_from_heme_to_add
#delete the now empty STC selection

alter (/6HR0/A/A/CYS`58),resi=str(int(resi)+3000)
alter (/6HR0/A/A/CYS`61),resi=str(int(resi)+3000)
alter (/6HR0/A/A/HIS`62),resi=str(int(resi)+3000)
#change the heme residues so they are not overlapping down the line


z = 0
#set z to be the same as number of hemes
number_of_hemes = int(number_of_hemes)
#set y to be number_of_hemes

python
for x in range(number_of_hemes):
    cmd.align("6HR0", "cxxch_motifs")      
     #puts 6HR0 onto the protein by aligning
    cmd.set_name("6HR0", "protein_heme")
     #renames 6HR0 as protein_heme
    cmd.copy("6HR0", "protein_heme")
     #copies 6HR0 to make a new heme for aligning
    cmd.delete("cxxch_motifs")
     #removes CXXCH motifs
    cmd.set_name("protein_heme", "protein_heme" + str(z))
     #renames protein_heme with a number based on it's insertion rank protein_heme
    cmd.select("already_inserted", "all w. 5 of (elem Fe)") 
    cmd.select("cxxch_motifs", "((ch_motifs (or (r. CYS and (ch_motifs xt. 20)))) and not (already_inserted xt. 8))")
     #remakes CXXCH motifs selection but now not including anything 1.5 angstroms or closer to an already inserted heme
    z = z + 1
     #ups z value
    cmd.delete("already_inserted")
python end

delete 6HR0
#remove remaining left over heme from file

cmd.select("old_residues", "i. 3058 or (i. 3061 or (i. 3062))")
cmd.remove("old_residues")
cmd.delete("old_residues")
#gets rid of residues from STC left behind

cmd.select("chain_residues", ("c. %s"%G))
select counting_temp, chain_residues and name ca
delete chain_residues
b = cmd.count_atoms("counting_temp")
number_of_residues = b
print number_of_residues
delete counting_temp

j = 0
q = 0
w = b + 1
n = 0
#sets j and q to 0 for upcoming python "for" loops


python

for x in range(number_of_hemes):
    cmd.alter("protein_heme" + str(q), 'resi=str(w)')
    cmd.alter("protein_heme" + str(q), 'chain=str(G)')
     #alters the residue of all the hemes to make them follow on from the final one of the main chain
    w = w + 1
    q = q + 1

cmd.select("new_hemes", "all and hetatm")
cmd.create("heme_chain", "new_hemes")
cmd.remove("new_hemes")
cmd.delete("new_hemes")
#creates a selection out of the hemes called heme chain

for x in range(number_of_hemes):
    cmd.delete("protein_heme" + str(j))
     #removes a protein heme objectect
    j = j + 1
#deletes the old heme chain scaffold
python end

remove  r. CYS and n. HG and (all w. 5 of (all and HETATM))
#deletes the cysteine hydrogen of the cxxch motifs

cmd.delete("cxxch_motifs")
cmd.delete("ch_motifs")
cmd.delete("cxxch_histidines")
#clears up used selections

cmd.create("output", "all")
remove not output
delete not output

select sele123, hetatm
alter "sele123",chain=str(G)
delete sele123
#renames hemes to be of the specified chain

save output_hemes_attached.pdb, all