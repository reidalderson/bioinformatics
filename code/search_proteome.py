class Proteome(object):
    """ Search for motifs within proteomes and compare the observed vs. expected motif frequencies."""
    def __init__(self):
        """ 
        Store the amino acids in single-letter format for later use 
        """
        self.amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        
        
    def load_uniprot(self, uniprot):
        """ 
        Read a UniProt proteome file.
        
        Note that the input proteome file should contain the following column order:
        (1) Entry (2) Entry name (3) Status (4) Protein names (5) Gene names (6) Organism (7) Length (8) Sequence 
            
        Parameters
        ----------		
        `uniprot` : file downloaded from UniProt containing a proteome of interest
        
        Returns
        -------
        `uniprot_proteome` : pandas.DataFrame, contains the FASTA IDs, sequences, and length of sequences 
        """        
        # Load the dataset with column headers defined
        self.uniprot_proteome = pd.read_csv(uniprot, sep='\t', names=['Entry', 'Entry name', 'Status', 'Protein names', 'Gene names', 'Organism', 'Length', 'Sequence'])
        
        # Convert protein length values to integers
        self.uniprot_proteome['Length'][1:] =  self.uniprot_proteome['Length'][1:].astype(int)  
        
        return self.uniprot_proteome
        
 
    def load_fasta(self, fasta):
        """ 
        Read a FASTA file with identifiers and sequences.
        
        Note that the FASTA file can be single- or multi-line.
            
        Parameters
        ----------		
        `fasta` : file, standard FASTA format with identifiers (e.g. >) and sequences 
        
        Returns
        -------
        `fasta_proteome` : pandas.DataFrame, contains the FASTA IDs, sequences, and length of sequences
        """
        # Load the fasta file
        self.seq_list = SeqIO.index(fasta, "fasta")  # makes dictionary that contains Seq objects

        # Create arrays to store ID and sequence information
        ids = []  ; seqs = []
        for key, value in self.seq_list.items():
            ids.append(key)
            seqs.append(str(value.seq)) # convert the Seq object to a string

        # Create a dictionary for subsequent conversion into a dataframe
        self.dictionary = {'ID':ids, 'Sequence':seqs} 

        # Load the dictionary into a dataframe
        self.fasta_proteome = pd.DataFrame(self.dictionary)  
        self.fasta_proteome['Length'] = self.fasta_proteome['Sequence'].str.len()  # create column in the dataframe with sequence lengths
        self.fasta_proteome['Length'][1:].astype(int)   # convert the length values to integers

        return self.fasta_proteome
        
        
    def proteome_length(self, query_proteome):
        """ 
        Calculate the total number of amino acids in a proteome.
            
        Parameters
        ----------		
        `query_proteome` : pandas.DataFrame, contains the proteome of interest
        
        Returns
        -------
        `total_length` :  int, the total number of residues in the queried proteome
        """   
        # Sum up the total number of residues
        self.total_length = query_proteome['Length'][1:].sum(axis=0, skipna=True)
        
        return self.total_length
        
        
    def proteome_AA_fractions(self, query_proteome):
        """ 
        Calculate the fractional amino acid composition of a proteome.
            
        Parameters
        ----------		
        `query_proteome` : pandas.DataFrame, contains the proteome of interest
        
        Returns
        -------
        `AA_fractions` :  pandas.DataFrame, the amino acid frequencies in the proteome
        """   
        
        # Initialize an empty array to store the fractional amino acid values
        fractions = []
    
        # Compute the fractional amino acid composition of the proteome
        for aa in self.amino_acids:
            fractions.append(self.find_motifs(query_proteome, aa)[1])
            
        # Convert into a dictionary and then a Pandas dataframe 
        self.AA_fractions = pd.DataFrame(list(dict(zip(self.amino_acids, fractions)).items()), columns=['AA', 'Fraction'])
        self.AA_fractions['Fraction'][1:].astype(float) # convert the fractions to floats
        
        return self.AA_fractions  
 
 
    def find_motifs(self, data_frame, query):
        """ 
        Search a proteome for a specific amino acid or motif.
            
        Parameters
        ----------		
        `data_frame` : Pandas dataframe, contains the proteome of interest
        
        `query` :  str, the amino acid or motif to be searched
        
        Returns
        -------
        `sum` :  int, total number of hits for the queried residue/motif
        
        `fraction`:  float, sum divided by the total number of residues or motifs of the same size
        """  
        # Extract sequences from the pandas dataframe
        seqs = data_frame['Sequence']
                
        # Calculate the total number of residues in the queried proteome
        total_residues = self.proteome_length(data_frame)
        
        # Initialize an array
        hits_array = []
        
               
        # Loop over each sequence and count the number of times the queried motif shows up
        for protein in seqs:
            hits_array.append(protein.count(query))

        # total number of hits for the queried residue/motif, check if any hits are found
        self.sum = np.sum(np.asarray(hits_array))  
        if self.sum < 1:
            self.fraction = 0
        else:
            self.fraction = self.sum/(total_residues - len(query)+1)  # total number of hits divided by the total number of residues or motifs of the same size
                        
        return self.sum, self.fraction       
 
        
    def find_flex_motifs(self, data_frame, list_of_motifs):
        """ 
        Search a proteome for a list of motifs.
            
        Parameters
        ----------		
        `data_frame` : int, contains the proteome of interest
        
        `list_of_motifs` :  list of str, the amino acids or motifs to be searched.
        
        Returns
        -------
        `counted_motifs` :  pandas.DataFrame, the queried motifs and number of instances of each motif.
        """         
        # Initialize arrays
        query_array = []  ; found_array = []  ; length_array = []
        
        # Loop over the list of motifs and count the number of instances
        for i, motif in enumerate(list_of_motifs):
            query_array.append(motif)
            found_array.append(self.find_motifs(data_frame, motif)[0])
            length_array.append(len(motif))
        
        # Create an array with the total number of possible motifs for each motif
        all_motifs = np.zeros(shape=(len(length_array),))
        for x in range(0, len(length_array)):
            all_motifs[x] = self.proteome_length(data_frame) - length_array[x]

        # Convert into a dictionary and then a Pandas dataframe 
        self.counted_motifs = pd.DataFrame(list(dict(zip(query_array, found_array)).items()), columns=['Motif', 'Count'])
        self.counted_motifs['Count'][1:].astype(int)  # make the Count column integer format         
        self.counted_motifs['Fraction'] = np.asarray(found_array)/np.asarray(all_motifs)  # add a column Fraction for the observed # motifs divided by the total # motifs
    
        return self.counted_motifs

        
    def iterate_motif(self, query):
        """ 
        Generate a list of all possible motifs that are consistent with an input motif 
        
        Note the following usage specifications:  the "_" symbol separates positions in the motif; "X" designates any residue;        
        and multiple residues within a given motif position indicate "or".
        
        Example:  MILV_X_DE means Met, Ile, Leu, or Val followed by any residue followed by Asp or Glu.
        
        Parameters
        ----------		
        `query` : str, the amino acid motif to be queried 
                use "_" to separate positions in the motif, e.g. I/V-X-I/V would be IV_X_IV
        
        Returns
        -------
        `all_possibilities` :  list of str, all variations of the inputted motif (e.g. IV_X_IV yields [IAI, ICI, IDI, ... , VVV, VWV, VYV]
        """ 
        # Split the input query into its corresponding residue positions
        listq = [item for item in query.split("_")]
        
        
        # If no "_" is supplied (i.e. IXI), then create list of residues
        if len(listq) == 1:
            listq = list(listq[0])
        else:
            pass

        # Generator an iterable list of motifs  
        positions = [[] for q in listq]
        for i in range(len(listq)):
            if len(listq[i]) == 1:
                if listq[i] == "X":
                    for aa in self.amino_acids:
                        positions[i].append(aa)
                else:
                    positions[i].append(listq[i])
            elif len(listq[i]) > 1:
                for j in listq[i]:
                    positions[i].append(j)
            else:
                print("WARNING! iterate_motif does not recognize your input")
                print('Please use single-letter amino acid codes (A, C, D, ...) or X to denote any amino acid')
                print('Note that residue positions should be separated by the "_" character (no " ")')
                print('exiting now....')
                sys.exit()
        
        # Use itertools to generate the final list 
        self.all_possibilities = [''.join(x) for x in list(itertools.product(*positions))]

        return self.all_possibilities
        
    def expected_motifs(self, query_proteome, query_motif):
        """ 
        Calculate the expected number of a queried motif based on amino acid frequencies in the proteome.
            
        Parameters
        ----------		
        `query_proteome` : Pandas dataframe
            Contains the proteome of interest.
        
        `query_motif` : str or list
            Contains the motif(s) of interest.
        
        Returns
        -------
        `final_dataframe` :  pandas.DataFrame, expected fraction and count of a given motif. 
        """   
        # Get fractional values of the amino acids within the proteome
        fractional_AA = self.proteome_AA_fractions(query_proteome)
        
        # Get length of the queried proteome 
        length = self.proteome_length(query_proteome)
        
        # Create a dictionary with amino acid-to-number mapping
        numbers = {a: b for b, a in dict(enumerate(self.amino_acids)).items()}
        
        # Check if query_motif is a list of strings
        if isinstance(query_motif, list) == True:
            probabilities = np.zeros(shape=(len(query_motif)))  # Create a numpy array to store probabilities for each motif
            expected = np.zeros(shape=(len(query_motif)))  # Create a numpy array to store expected total for each motif
            observed = np.zeros(shape=(len(query_motif)))  # Create a numpy array to store observed total for each motif
            for i, motif in enumerate(query_motif):
                for res_num in range(0, len(motif)):
                    if res_num == 0:
                        probabilities[i] = fractional_AA.at[numbers[motif[res_num]], 'Fraction']  # get the value at the specified position
                    else:
                        probabilities[i] *= fractional_AA.at[numbers[motif[res_num]], 'Fraction'] # multiply by the value at the specified position
                expected[i] = int(probabilities[i] * (length - len(motif)-1))  # multiply the probability and the total number of motifs
                observed[i] = int(self.find_motifs(query_proteome, motif)[0])
                
            # Combine into final Pandas dataframe    
            self.final_dataframe = pd.DataFrame({'Motif':query_motif, 'Probability':probabilities, 'Observed':observed, 'Expected':expected})
        
        
        # Check if query_motif is a string
        elif isinstance(query_motif, str) == True:
            probabilities = np.zeros(shape=(1,))  # Create a numpy array to store probabilities for each motif
            expected = np.zeros(shape=(1,))  # Create a numpy array to store expected total for each motif
            observed = np.zeros(shape=(1,))  # Create a numpy array to store observed total for each motif
            for res_num in range(0, len(query_motif)):
                if res_num == 0:
                    probabilities[0] = fractional_AA.at[numbers[query_motif[res_num]], 'Fraction']  # get the value at the specified position
                else:
                    probabilities[0] *= fractional_AA.at[numbers[query_motif[res_num]], 'Fraction'] # multiply by the value at the specified position
                expected[0] = int(probabilities[0] * (length - len(query_motif)-1))  # multiply the probability and the total number of motifs
                observed[0] = int(self.find_motifs(query_proteome, query_motif)[0])
                
            # Combine into final Pandas dataframe    
            self.final_dataframe = pd.DataFrame({'Motif':query_motif, 'Probability':probabilities, 'Observed':observed, 'Expected':expected})
            
        self.final_dataframe['Observed'][1:].astype(int)  # make the observed column integer format 
        self.final_dataframe['Expected'][1:].astype(int)  # make the expected column integer format
        
        return self.final_dataframe
        
        
    def find_proteins(self, data_frame, query, proteome_type):
        """ 
        Search for proteins in a proteome that contain the queried motif.
            
        Parameters
        ----------		
        `data_frame` : Pandas dataframe, contains the proteome of interest
        
        `query` :  str or list of strings, the amino acid or motifs to be searched (e.g. 'IPV' or ['IPI', 'IPV', 'VPI', 'VPV'])
        
        `proteome_type` : str, either 'FASTA' or 'UniProt' depending on the analyzed proteome
        
        Returns
        -------
        `found_proteins` :  list, array of UniProt/FASTA IDs that contain the queried residue/motif
        """  
        # First check which type of proteome/Pandas dataframe is being queried
        if proteome_type.upper() == 'FASTA':
            id = 'ID'
        elif proteome_type.upper() == 'UNIPROT':
            id = 'Entry'
        else:
            print('ERROR! The proteome_type %s is not recognized in the function "find_proteins" !' % proteome_type)
            print('Please enter "fasta" or "uniprot", or modify the function to accommodate other types')
            print('Exiting now.........................................................................')
            sys.exit()

        # Initialize an empty array
        self.found_proteins = []
        
        # Check if queried motif (i.e. query) is a list
        if isinstance(query, list) == True:
            for i in range(0, len(data_frame['Sequence'])):
                for j, motif in enumerate(query):
                    if motif in data_frame['Sequence'][i]:
                        self.found_proteins.append(data_frame[id][i])
              
        # If not list, check if query is a string
        elif isinstance(query, str) == True:
            for i in range(0, len(data_frame['Sequence'])):
                if query in data_frame['Sequence'][i]:
                    self.found_proteins.append(data_frame[id][i])

        else:
            print('In the function "find_proteins", the "query" must be either a list of strings or a single string....!')
            print('Please check what you entered!')
            print('Exiting now...................')
            sys.exit()
        
        # Remove duplicates from list
        self.found_proteins = list(dict.fromkeys(self.found_proteins))
        
        return self.found_proteins
        
        
    def make_motifs(self, query):
        """ Return motifs consistent with flex input"""
        self.comb_list = []
        for i in range(len(query.split("_"))):
            self.comb_list.append([])
            for j in query.split("_")[i]:
                self.comb_list[i].append(j)

        all_motifs = list(itertools.product(*self.comb_list))
        self.motifs = [''.join(x) for x in all_motifs]
        
        return self.motifs
       
       
    def calc_stats(self, query_proteome, query, difference='n', sub_proteome=None):
        """ """
        fname = '%s_observed_expected.xlsx' % str(query)

        motifs = self.make_motifs(query)
        columns = ['Obs_', 'Exp_']
        headers = ([j+i for i in motifs for j in columns])
        headers.insert(0, 'AA')
                
        # Prepare list to store data for the columns of a new DataFrame
        data = [self.amino_acids]       
        
        # Iterate over each motif and store the Observed/Expected frequencies
        for motif in motifs:
            if difference.upper() == 'N':
                calc = self.expected_motifs(query_proteome, self.iterate_motif(motif))
                data.append(calc['Observed'])
                data.append(calc['Expected'])
            else:
                try:
                    calc_main = self.expected_motifs(query_proteome, self.iterate_motif(motif))
                    calc_sub = self.expected_motifs(sub_proteome, self.iterate_motif(motif))
                    data.append(calc_main['Observed'] - calc_sub['Observed'])
                    data.append(calc_main['Expected'] - calc_sub['Expected'])
                except ValueError:
                    print('ERROR! Please include a 2nd proteome (sub_proteome) in the function "calc_stats" ')
                    print('Exiting now.........................................................................')
                    sys.exit()
        
        # Create the DataFrame using the newly created lists of Headers and Data
        final = pd.DataFrame.from_dict(dict(zip(headers, data)))
        
        # Calculate chi-squared values 
        n_motifs = len(motifs)
        df_rows = len(final.index)
        self.chisq_array = np.zeros(shape=(df_rows, n_motifs))

        # loop over AA types, over motifs
        for i in range(0, df_rows):
            for j in range(0, n_motifs):
                if j == 0:
                    exp = final.values[:,j+1]
                    obs = final.values[:,j+2]
                else:
                    exp = final.values[:,2*j+1]
                    obs = final.values[:,2*j+2]
                self.chisq_array[i,j] = ((obs[i]-exp[i])**2)/exp[i]

        # Calculate P values
        motif_chisq = np.sum(self.chisq_array, axis=0)                   # sum the chi-squared values down each column (i.e. for each motif)
        self.aa_chisq = np.sum(self.chisq_array, axis=1)                 # sum the chi-squared values across each row (i.e. for each amino acid)
        self.df = n_motifs*df_rows#len(motifs[0])                        # degrees of freedom, i.e. the number of residues in the motif


    def plot_chisq(self, query, pval=0.01, diff='n', main=None, sub=None, title=None, save_flg='y'):

        if diff.upper() == 'N':
            bar_vals = self.aa_chisq
        else:
            bar_vals = self.calc_stats(main, query, diff='y', sub_proteome=sub)
            bar_vals = self.aa_chisq
            
        plt_array = []
            
        # make summed chisq bar graph
        fig = plt.figure()
        bar = plt.bar(self.amino_acids, bar_vals, color='k')
        plt.ylabel(r'$\sum_{i=1}^%.0f \chi^{2}_{i}$ value' % len(self.make_motifs(query)))
        xlab = []
        for i, res in enumerate(query.split("_")):
            if len(res)>1:
                new = '['
                new += '/'.join(list(res))
                new += ']'
                xlab.append(new)
            else:
                xlab.append(res)
        xlab = '-'.join(xlab)
        plt.xlabel('Amino acid in %s' % xlab)
        if not title:
            plt.title(r'$\chi^{2}$ values for %s' % (query))
        else:
            plt.title(r'$\chi^{2}$ values for %s in %s' % (query, title))
         
        df = 0
        for res in query.split("_"):
            if len(res)>1:
                df += len(list(res))
            else:
                df += 1
        df = df-1
        print('>> There are %.0f degrees of freedom' % df)
        
        
        # make stacked bar graph
        fig = plt.figure()
        labels = self.make_motifs(query)
        colors = ['k', 'r', 'purple', 'blue', 'darkgreen', 'orange', 'brown']*2
        for i in range(np.shape(self.chisq_array)[1]):
            if i < np.shape(self.chisq_array)[1]-1:
                if i == 0:
                    plt.bar(self.amino_acids, self.chisq_array[:,i], label=labels[i], color=colors[i])
                else:
                    plt.bar(self.amino_acids, self.chisq_array[:,i], bottom=sum(self.chisq_array[:,j] for j in range(0,i)), label=labels[i], color=colors[i])
            else:
                bar = plt.bar(self.amino_acids, self.chisq_array[:,i], bottom=sum(self.chisq_array[:,j] for j in range(0,i)), label=labels[i], color=colors[i])

        # Add counts above the two bar graphs
        for i, rect in enumerate(bar):
            height = rect.get_height() + np.sum(self.chisq_array[i,:-1])
            if stats.chi2.sf(height, df) < pval: 
                print('>>> yes! ', i, stats.chi2.sf(height, df))
                plt.text(rect.get_x() + rect.get_width()/2.0, height, '*', ha='center', va='bottom') 
                
        plt.legend(loc='upper right')
        plt.ylabel(r'$\chi^{2}$ value')
        plt.xlabel('Amino acid in %s' % xlab)
        if not title:
            pass
        else:
            plt.title('%s' % title)
        plt.tight_layout()    
        if save_flg.upper() != 'Y':
            plt.show()  
        else:
            plt_array.append(fig)
        
        return plt_array
