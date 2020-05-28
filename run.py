from search_proteome import Proteome 
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

#########################################################################################
##########      Usage example of the search_proteome.py script                 ##########
##########            Imports the Proteome() class                             ##########
###########                  Reid Alderson, May 2020                           ##########
#########################################################################################


if __name__ == "__main__":
    print('\n**** Beginning example demo ****\n')

    #### Initialize the Proteome Class for testing      
    start = Proteome()

    #### Example (1) Load the proteomes with the functions "load_uniprot" and "load_fasta", depending on input file type       
    uniprot = start.load_uniprot('input/uniprot_human.txt')
    idrs = start.load_fasta("input/human_idrs.fasta")
    print('###############################################################################')
    print('\n(1) Loaded proteome and IDR-ome successfully!\n')
    print('###############################################################################\n')


    ### Example (2) Compute the length of the two proteomes with the function "proteome_length"
    print('###############################################################################')
    print('\n(2) Calculating total length of the proteomes')
    print('\nLength of proteome = ', start.proteome_length(uniprot))  
    print('Length of "IDR-ome" = ', start.proteome_length(idrs))
    np.savetxt('output/2_proteome_length.txt', np.asarray([start.proteome_length(uniprot)]), fmt='%s', header='Total number of residues in the proteome')
    np.savetxt('output/2_disordered_length.txt', np.asarray([start.proteome_length(idrs)]), fmt='%s', header='Total number of residues in the IDR-ome') 
    print('\n###############################################################################\n')


    ### Example (3) Calculate amino acid frequencies in the proteomes with the function "proteome_amino_acid_fractions"
    print('###############################################################################')
    print('\n(3) Calculate amino acid frequencies in the proteome')
    print('\nAA fractions in the proteome ---->\n', start.proteome_amino_acid_fractions(uniprot))
    print('\nAA fractions in the "IDR-ome" ---->\n', start.proteome_amino_acid_fractions(idrs))
    np.savetxt('output/3_disordered_AA_frequencies.txt', start.proteome_amino_acid_fractions(idrs).values, fmt='%s', header='Amino acid frequencies in the IDR-ome')
    np.savetxt('output/3_proteome_AA_frequencies.txt', start.proteome_amino_acid_fractions(uniprot).values, fmt='%s', header='Amino acid frequences in the human proteome')
    print('\n###############################################################################\n')



    ### Example (4) Count the number of specific motifs with the function "find_motifs"
    query = 'IPV'
    print('###############################################################################')
    print('\n(4) Counting specific motifs in the proteome')
    print('\n%s motifs in proteome = %s' % (query, start.find_motifs(uniprot, query)))
    print('%s motifs in "IDR-ome" = %s' % (query, start.find_motifs(idrs, query)))
    np.savetxt(('output/4_disordered_%s_motif_count.txt' % query), start.find_motifs(idrs, query))
    np.savetxt(('output/4_proteome_%s_motif_count.txt' % query), start.find_motifs(uniprot, query))
    print('\n###############################################################################\n')



    ### Example (5) Count the number of a list of motifs with the function "count_multiple_motifs"
    list_motifs = ['IPI', 'IPV', 'VPI', 'VPV']
    print('###############################################################################')
    print('\n(5) Counting multiple motifs in the proteome')
    print('\n%s motifs in proteome =\n%s' % (list_motifs, start.count_multiple_motifs(uniprot, list_motifs)))
    print('\n%s motifs in "IDR-ome" =\n%s' % (list_motifs, start.count_multiple_motifs(idrs, list_motifs)))
    np.savetxt(('output/5_disordered_%s_motifs_count.txt' % '-'.join(list_motifs)), start.count_multiple_motifs(idrs, list_motifs).values, fmt='%s')
    np.savetxt(('output/5_proteome_%s_motifs_count.txt' % '-'.join(list_motifs)), start.count_multiple_motifs(uniprot, list_motifs).values, fmt='%s')
    print('\n###############################################################################\n')



    ### Example (6) Count the number of motifs using a variable input with the functions "iterate_motif" and "count_multiple_motifs"
    query = 'IV_X_IV'
    iterator = start.iterate_motif(query)
    print('###############################################################################')
    print('\n(6) Counting multiple motifs in the proteome using "iterate_motif"\n')
    print('Proteome = \n', start.count_multiple_motifs(uniprot, iterator))
    print('"IDR-ome" = \n', start.count_multiple_motifs(idrs, iterator))
    np.savetxt(('output/6_disordered_%s_motifs_count.txt' % query), start.count_multiple_motifs(uniprot, iterator).values, fmt='%s')
    np.savetxt(('output/6_proteome_%s_motifs_count.txt' % query), start.count_multiple_motifs(idrs, iterator).values, fmt='%s')
    print('\n###############################################################################\n')


    ### Example (7a) Now write out the expected number of motifs based on AA frequencies compared to the observed total 
    query = 'IV_X_IV'
    iterator = start.iterate_motif('I_X_V')
    print('###############################################################################')
    print('\n(7) Save the expected number of motifs based on AA frequencies compared to the observed total with "expected_motifs"\n')

    print('Proteome = \n', start.count_multiple_motifs(uniprot, iterator))
    print('"IDR-ome" = \n', start.count_multiple_motifs(idrs, iterator))
    np.savetxt(('output/7a_proteome_%s_expected_observed.txt' % query), start.expected_motifs(uniprot, start.iterate_motif(query)).values, fmt='%s')
    np.savetxt(('output/7a_disordered_%s_expected_observed.txt' % query), start.expected_motifs(idrs, start.iterate_motif(query)).values, fmt='%s')
    print('\n###############################################################################\n')


    ### Example (7b) Now make plots

    # first set up matplotlib parameters
    import matplotlib as mpl
    mpl.rcParams['axes.linewidth'] = 0.5           # linewidth on axes
    mpl.rcParams['axes.titlesize'] = 14            # fontsize of title
    font_style = 'Arial'                           # fontstyle
    mpl.rcParams['font.family']=font_style         # set font
    mpl.rcParams['mathtext.fontset'] = 'stix'      # stix font is readable by Adobe Illustrator version CC
    mpl.rcParams['pdf.fonttype'] = 42              # use this to output PDF with modifiable text
    mpl.rcParams['xtick.major.size'] = 5           # x-tick length
    mpl.rcParams['xtick.major.width'] = 0.5        # x-tick linewidth
    mpl.rcParams['ytick.major.size'] = 5           # y-tick length
    mpl.rcParams['ytick.major.width'] = 0.5        # y-tick linewidth    
    mpl.rc('axes', labelsize=14)                   # fontsize of axis labels

    mpl.rc('xtick', labelsize=12)                  # fontsize of the xtick labels
    mpl.rc('ytick', labelsize=12)                  # fontsize of the ytick labels

    # now make the plots
    start.calc_stats(uniprot, query, difference='y', sub_proteome=idrs)
    structured_plot = start.plot_chisq(query, pval_cutoff=0.01, difference='y', main_proteome=uniprot, sub_proteome=idrs, plt_title='Structured', save_flg='y')
    plt.savefig('output/7b_structured_chisq.pdf')
    plt.savefig('output/7b_structured_chisq.png')
    start.calc_stats(idrs, query, difference='n')
    disordered_plots = start.plot_chisq(query, pval_cutoff=0.01, difference='n', main_proteome=idrs, plt_title='Disordered')
    plt.savefig('output/7b_disordered_chisq.pdf')
    plt.savefig('output/7b_disordered_chisq.png')

        
        
    ### Example (8) Find proteins that contain a certain motif using the function "find_proteins"
    print('###############################################################################')
    print('\n(8) Find proteins that contain a certain motif using the function "find_proteins"')
    query = 'IPV'
    np.savetxt(('output/8_disordered_protein_%s_motifs.txt' % query), start.find_proteins(idrs, query, 'fasta'), fmt='%s', header=('ENSEMBL IDs of proteins in the IDR-ome that contain at least one %s motif' % query))
    np.savetxt(('output/8_proteome_proteins_%s_motifs.txt' % query), start.find_proteins(uniprot, query, 'uniprot'), fmt='%s', header=('UniProt IDs of proteins in the proteome that contain at least one %s motif' % query))
    print('\n###############################################################################\n')


    ### Example (9) Now calculate the total number of [I/V]-X-[I/V] motifs in structured vs. disordered proteomes 
    print('###############################################################################')
    print('\n(9) Now calculate the total number of [I/V]-X-[I/V] motifs in structured vs. disordered proteomes \n')
    query = 'IV_X_IV'                                       # load the query motif
    iterator = start.iterate_motif(query)                   # get all motifs 
    data_array = np.zeros(shape=(2, len(iterator), 2))      # set up 3D data array 

    # store the data in a numpy array and write to disk
    fout1 = open(('output/9_proteome_count_%s_motifs.txt' % query), 'w')
    fout2 = open(('output/9_disordered_count_%s_motifs.txt' % query), 'w')
    fout1.write('%s\t%s\t%s\n' % ('Motif', 'Count', 'Fraction'))
    fout2.write('%s\t%s\t%s\n' % ('Motif', 'Count', 'Fraction'))

    # first store the data, calculate the total count, and write out the total
    for i, motif in enumerate(iterator):
        data_array[0, i, 0], data_array[0, i, 1] = start.find_motifs(uniprot, motif)
        data_array[1, i, 0], data_array[1, i, 1] = start.find_motifs(idrs, motif)  
    fout1.write('%s\t%.0f\t%f\n' % ('Total', np.sum(data_array[0,:,0]), np.sum(data_array[0,:,1])))
    fout2.write('%s\t%.0f\t%f\n' % ('Total', np.sum(data_array[1,:,0]), np.sum(data_array[1,:,1])))

    # then write out the count for each motif 
    for j, val in enumerate(iterator):
        fout1.write('%s\t%.0f\t%f\n' % (val, data_array[0,j,0], data_array[0,j,1]))
        fout2.write('%s\t%.0f\t%f\n' % (val, data_array[1,j,0], data_array[1,j,1]))
    fout1.close()
    fout2.close()


    # set up plot 
    fig, ax = plt.subplots(2, 4, sharey=True, sharex=True, figsize=(11,4))
    ax = ax.ravel() 

    # create list of motifs to iterate over and specify the location of subplots within the grid 
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    motif_list = ['I_X_I', 'I_X_V', 'V_X_I', 'V_X_V']
    structured_list = [0, 1, 4, 5]
    structured_dict = dict(zip(motif_list, structured_list))
    disordered_dict = dict(zip(motif_list, [x+2 for x in structured_list]))

    # Iterate over the list of motifs and calculate proteome-disordered = structured, plot structured vs. disordered 
    for i, motif in enumerate(motif_list):
        proteome_data = start.count_multiple_motifs(uniprot, start.iterate_motif(motif)).to_numpy()
        idr_data = start.count_multiple_motifs(idrs, start.iterate_motif(motif)).to_numpy()
        ax[structured_dict[motif]].bar(amino_acids, (proteome_data[:,1]-idr_data[:,1])/(start.proteome_length(uniprot) - start.proteome_length(idrs) - len(motif)+1), color='k' )
        ax[disordered_dict[motif]].bar(amino_acids, idr_data[:,2], color='r')
        
    # Set axes labels & titles
    ax[0].set_ylabel('% of tripeptides')
    ax[4].set_ylabel('% of tripeptides')
    ax[4].set_xlabel('X residue type')
    ax[5].set_xlabel('X residue type')
    ax[6].set_xlabel('X residue type')
    ax[7].set_xlabel('X residue type')
    ax[0].text(0.8, 1.1, 'Structured', transform=ax[0].transAxes, size=14, va='center', ha='center')  
    ax[2].text(0.8, 1.1, 'Disordered', transform=ax[2].transAxes, size=14, va='center', ha='center')

    # Label each subplot
    labels = ['IxI', 'IxV', 'IxI', 'IxV', 'VxI', 'VxV', 'VxI', 'VxV']
    for n, axs in enumerate(ax):   
        axs.text(0.02, 0.88, labels[n], transform=axs.transAxes, size=12)

    # Fix the layout
    plt.tight_layout(h_pad=-0.05, w_pad=-0.7)
    plt.savefig('output/9_compare_fractions_%s_motifs.pdf' % query)

    print('###############################################################################\n')   
    print('\n**** Finished example demo ****\n')
    print('\n**** Please find the output files in the output directory ****\n\n')    
