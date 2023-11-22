# PROTEIN-PROTEIN INTERFACE ANALYSIS 
Receptor Binding Domain (RBD) of SARS-Cov-2 Spike Protein and its receptor Angiotensin Converting Enzyme (ACE2)


Carmen-Juliette Samedi
Gemma Sagastegui
Júlia Prenafeta

INDEX


INTRODUCTION                                                                                                 3
Objectives
Background 
Strategy 


METHODOLOGY
Preparation step                                                                                             4 
Step 1                                                                                                       5
Step 2                                                                                                       8
Step 3                                                                                                       9
Step 4                                                                                                       10

DISCUSSION ON RESULTS 
Problems encountered                                                                                         11
Conclusion                                                                                                   12










Introduction

Objective
To appraise the relative contribution of interface residues to the interaction energy in a protein-protein complex. In our case, we are facing the interaction between the Receptor Binding Domain (RBD) of SARS-Cov-2 Spike Protein and its particular receptor, Angiotensin-converting enzyme 2.

Background

Proteins rely on interactions with other molecules for their function (often vital so as to control their activity).  It is important to understand how changes in amino acid sequences impact these interactions, by remarking the need of studying how amino acids at the protein-protein interface influence the energy involved in their interaction within a complex. 

In this specific case study, we will delve into the interaction between the Receptor Binding Domain (RBD) of the SARS-CoV-2 Spike protein and its receptor, ACE2. This interaction is crucial in viral infection. The analysis is centered on the 6m0j structure from the PDB, which only contains the RDB domain.

We will also use a significant metric in this study: ΔΔG (Delta Delta G), which predicts how a single point mutation affects protein stability by measuring the energy change between folded and unfolded states. The prediction informs us about the favorability of a variant concerning protein stability.

Strategy
Evaluation of the contribution of individual amino acids residues through the succeeding steps:
Get to know amino acid residues that form the interface between the complex components (by using the PDB molecule 6m0j)
Compute the stability contribution to the complex by following a traditional Ala-scanning experiment (for replacing each residue in turn by Ala and evaluating the changes in the complex interaction energy.
Identify known SARS-Cov-2 Spike Variants and evaluate their effect on ACE2 Binding, plus replacing the appropriate residue with the variant and reevaluate the interaction energy

Methods
The preparation previously starting is necessary so as to obtain the two clean molecules of interest, by means of without residues around them. 
We first downloaded the protein structure from PDB (6m0J) and removed all chains except the ones involved in the biological unit. 
Afterwards, we carried out a quality check on the structures. This implied adding missing hydrogen atoms, the structures and atom charges, but also, the removal of the number of heteroatoms, water molecules and ligands. Using the code provided on github by Josep-Lluis we were able to do the following,
Preparation step 

The structure before cleaning the protein, many Heteroatoms, many water molecules, ligands and metal residues but no hydrogen atoms.

After doing the quality check, all heteroatoms were removed, no more water molecules and almost 100 residues were removed in total. Hydrogens were added. 



















Step 1
The aim of this first step was the visualization of the 2 structures, its differentiation and also, the selection of the atoms involved in the interface and finding polar bonds. 
We used Pymol for all of the requirements mentioned above.

We actually found 2 ways of finding the same interactions between chains:
1st way: by directly finding interactions between chains.

Then, showing its distance












2nd way: by selecting first the possible areas with interactions between chains


Then, finding the polar contacts to any atoms.



After that, we showed both chains as lines so we can select the residues of the interactions we found easier, but only the ones between chains. Yellow lines shown are the interactions.










Then, we only showed as sticks the residues we have selected. Now, we are ready to use the wizard tool.


In order to use wizard, we have to click on the upper toolbar and then find ‘Measurements’. Then, we selected the atoms from each interaction to measure its distance.



Result obtained:












We obtained, and therefore, concluded that the largest distance is 2.6 Å (between LYS 31 and GLN 493), and so as to ensure that the adjacent residues are also considered, we added 1.4 Å, giving us the result of 4 Å as the largest distance.
We also prepare a python script (step1-interface-res.ipynb) to look for the residues that have contact between chains. This are the ones we obtained:

Contacts between residues in different chains:
GLN 42 A - TYR 449 E - 2.9429228
LYS 353 A - GLY 502 E - 2.7658923
TYR 41 A - THR 500 E - 2.624651
GLU 35 A - GLN 493 E - 2.7958221
TYR 83 A - ASN 487 E - 2.6658735
HIE 34 A - TYR 453 E - 2.8842824
ASP 30 A - LYS 417 E - 2.9862466
GLN 24 A - ASN 487 E - 2.615844
ASP 38 A - TYR 449 E - 2.6912184

Step 2
The posterior step after setting the largest distance was composed by the evaluation of the interaction energy between chains in the bound state (the complex) and the unbound state (the two chains isolated in solution). 
We previously assumed that 3D structure does not vary between bound and unbound states. Therefore, bonded/non-bonded terms between atoms from the same chain won’t change, and consequently, won’t be considered.
We also considered solvation energies from the ASA values for all atom types (not only the hydrophobic ones).

As a formula to compute the binding energy we used the difference of electrostatic, Van der Waals and solvation energy between chains, minus the solvation energy individually of each of them. 



By using Biopython, we went through the protein structure so as to gather chains and residues IDs. With the help of the variable called MAXDIST, which we set to 4 Å (as we previously mentioned), we computed various values by calling specific functions. Afterwards, the formula ensured us to use these values so as to find the total energy.
We performed this experiment with all the residues, and the interface ones,  too. 
These two situationships were different with the help of the variable MAXDIST, as in the all residues situation, it was set to zero. This guaranteed us that the formula computed the energy for ALL residues in the protein.

Final ΔG=-172.639878968562 for residues and ΔG=-127.88949147947226 for residues 
For all residues:				For the interface residues:
Step 3
Within the following step, we needed to determine the effect of replacing each interface residue by Ala in the overall ΔGA-B , and plot the results, highlighting the most relevant residues for the stability of the interface. We computed it by using the code of below:

# Effect of replacing each interface residue with Alanine
with open("res_ala.txt", "a") as res_ala:
   for ch in st[0]:
       for res in ch.get_residues():
           # Again, if cut off distance is bigger than 0 and residue is in the interface
           # calculate the following, if not continue
           if MAXDIST > 0 and res not in interface[ch.id]:
               continue
           print('{:1}, {:1.4}'.format(residue_id(res),
                   - elec[res] + elec_ala[res] - vdw[res] + vdw_ala[res] -solvAB[res] +\
                       solvAB_ala[res] -solvA[res] + solvA_ala[res]), file = res_ala)



Previously to running the previous code, we had to determine a variable with the possible atom names that could correspond to Ala atoms, which is the following:

#Possible Atom names that correspond to Ala atoms"
ala_atoms = {'N', 'H', 'CA', 'HA', 'C', 'O', 'CB', 'HB', 'HB1', 'HB2', 'HB3', 'HA1', 'HA2', 'HA3'}

We obtained the following results:


















Here you can observe the plot of the results obtained:



Step 4
Pymol allowed us to obtain images of the interface highlighting the relevant residues, and also, the interactions.
We chose the following colors so as to differentiate the different impacts values and the ones which had not.
The residues in color red : those with highest impact. Values between 5 and + 10 in DG difference.
The residues in color orange : those with the lowest impact. Values between 0.1- 4.9 in DG difference
The residues in color gray : those without impact. We later noticed that they were mostly composed of GLY (glycines).














Results and discussion
Problems
We encountered a problem once our programs were done; our results didn’t make sense. At first we reached out to Irene and realized it might have been an error in our distance measurements, however it came to be that we were using an erroneous pdb file of the protein. After looking back on the preparation step and what we did we realized we used a certain way to remove heteroatoms in Pymol which affected the result of the rest of the steps… 

OBTAINED                                           VS                            EXPECTED 

























Conclusions
Once we got all the results by following the steps, we obtained some conclusions.
To sum up, we have been able to determine how changes in the building blocks of a protein affect its stability, by using ΔG values to show it. 

We now focus on STEP 3. The obtained results show us the different ΔG of the change in the residue in the first column for an Alanine atom.  If we take as an example THR A27 (located in the second row of the first column), by changing it to an Ala atom we obtained a ΔG of 3.783. Alanine is a hydrophobic amino acid, and this property will help us to conclude some of the results.

By now taking into account the nature of the involved amino acids, we should previously mention that the Alanine atom side-chain is part of all others, but not for Glycine (GLY), and there won’t be need to do the replacement, and then, we had to get in mind only the different atoms.
We first start by looking at the distance in which the amino acid is from another's chain, as depending on the distance, if we replace the amino acid with an Alanine, it can be shorter, and then, they can not interact as they won’t reach the contact needed. Given this, some amino acid changes can either stabilize or destabilize the protein. 

If we take a closer look at the binding energy values, we can see that they can be differentiated between positive and negative values. For example, the ones obtained with the amino acids TYR (either A41, A83) or GLN (A42, E493), will destabilize the protein when changing for an Ala atom, as the binding energy is positive.
In one hand we have GLN, which is polar and then, hydrophilic, it will create some contacts with the H2O surrounding the structure, and it is exactly the contrary Ala will do, as previously mentioned, it is hydrophobic, so it won’t create contacts with H2O, destabilizing the protein.
On the other hand we have TYR, which gives us either positive (TYR A41, with a binding energy of 3.585) or negative (TYR E449, with a binding energy of -10.57). This can be explained with this amino acid’s structure. It has an aromatic ring, which makes it more hydrophobic, but also contains a hydroxyl group, which also makes amino acids more hydrophilic. 
Negative values obtained, such as the TYR E449 situation, as Ala is less hydrophobic than it, so Alanine will create contacts with H2,O, and will make the protein more stable.

We concluded that, when changing a polar amino acid (for example, Lysine (LYS) or Aspartic Acid (ASP)) by a nonpolar one, which is Alanine (ALA), it, will destabilize the protein, as a result of missing bonds with water that polar residues do.
Then, when this change is applied with a nonpolar amino acid (for example, PHE E486), it will stabilize the protein more. Phenylalanine is more non-polar than Alanine, as it side chain is a phenyl group, which is a large hydrophobic structure and then, doesn’t interact quite well with H2O, while Alanine has a simpler side chain, as it contains a methyl group (also hydrophobic but less in comparison to the phenyl group.

In other cases, we found some amino acids without or with little impact and change of the binding energy. Afterwards, we realized that it may be a consequence of distance, which can be either more or less equal, and then, it won’t affect the stability of the protein.

If we put all of our attention into STEP 4, the values we can obtain by swapping residues in the first column are either positive, negative or zero. We could reason that the more negative value we obtained, the more stable the protein became.

Afterwards, with the help of a graph, we could visualize which building blocks play the most significant role in maintaining the protein structure together.
In the y-axis will lay the binding energy and in the x-axis we put its ID. 
Given some residues were related with negative binding energy, we decided to put the binding energy in absolute values. This allowed us to see which are higher, and, consequently, those which had a higher impact on the structure of the protein. 
The performance of the alanine change allowed us to determine which residues are more relevant for the stability of the interface, and either the positive or negative effect on the protein.

We summarized the whole knowledge we obtained by analyzing this two complexes in 3 main points:
Structure. Pymol was crucial so as to see atoms' distance and how they connect with each other within the protein. We also have been aware that distances would differ, in other words, they won’t be exactly the actual distances, depending on the method of each examination of the structure they used.
Energy check. By comparing the total energy of all residues of the protein with the energy of the interface residues, and given both of them are negative values and quite similar, we could confirm our computations we have done are pretty accurate.
Low or high-impact blocks. Some block changes were not giving us a huge difference, as the distance between the blocks remained constant, the differences between the original and new blocks were not significant.


