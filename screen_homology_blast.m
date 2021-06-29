function [is_excluded_for_homology] = screen_homology_blast(query_sequence, refseq_database, gene_name, accession_number, params)
    % screen_homology_blast takes a query sequence (query_sequence), a path
    % to the RefSeq database, the current gene name and accession number of
    % interest, and the parameters struct. The function returns a boolean
    % variable excluded_for_blast that represents if the query sequence
    % passed the sequence homology specificity requirements.
    %
    % The relevant parameters used are: 
    % (1) e_value, which defines the E value used in the BLAST query.
    % (2) min_homology_length_to_exclude, which defines the minimum number
    % of matching bases in a BLAST hit.
    % (3) fivePrime_end_position and threePrime_start_position, used to
    % evaluate how well centered BLAST hits are.
    % (4) ligation_junction_overlap_cutoff, which defines the threshold of
    % overlapping bases on either side for centeredness.
    % 
    % The sequence specificity requirements are that: there are no BLAST
    % hits in other transcripts with the specified E value (or higher),
    % with more matching bases in the hit greater than or equal to
    % min_homology_length_to_exclude, that are also well-centered, i.e. 
    % span both sides of the ligation junction with at least
    % ligation_junction_overlap_cutoff bases on either side.
    
    e_value = params.E_VALUE;
    fivePrime_end_position = params.FIVEPRIME_END;
    threePrime_start_position = params.THREEPRIME_START;
    min_homology_length_to_exclude = params.BLAST_HOMOLOGY_MIN_LENGTH_CUTOFF;
    ligation_junction_overlap_cutoff = params.BLAST_HOMOLOGY_LIGATION_JUNCTION_MIN_OVERLAP_CUTOFF;
        
    is_excluded_for_homology = false;

    % Delete a fasta file left over from a premature termination of this
    % loop. If the program is run without deleting temp.fa, the BLAST
    % queries are not correct.
    delete('temp.fa');            
    fastawrite('temp.fa', query_sequence); 
    
    blastlocal_parameters = ['-S 1 -e ', num2str(e_value)];
    
    % Run BLAST query and save output in blastout. This assumes that the
    % BLAST has been added to the PATH variable.
    blastout = blastlocal('InputQuery','temp.fa','Program','blastn','DataBase', refseq_database, 'blastargs', blastlocal_parameters); 

    % Loop through all BLAST hits, examining if any of them are sufficient
    % to exclude the query sequence
    for ii = 1:length(blastout.Hits)              
       
        % The first step of the search is to ignore any hits in the same
        % gene, found by searching for the gene symbol and accession number
        % of the current gene of interest in the hit. If the gene symbol or
        % accession number are found, it is (by default) accepted as a
        % homology sequence.
        %        
        % If the gene symbol and accession number are not found, a careful
        % exclusion criteria search is performed.         
        %
        % Note that some genes have non-coding variant RNAs in RefSeq that
        % should also be ignored as hits. These transcripts typically start
        % with the gene symbol, and append a P and a number, representing a
        % pseudogene. The precise search query used to find matching gene
        % symbols is '(gene_name'. The leading parenthesis ensures that the
        % gene symbol of the gene of interest is at the start of the gene
        % symbol in the hit. For example, the search for KRT17 will result
        % in hits in KRT17P1, a non-coding RNA pseudogene for KRT17 (among
        % other KRT17 pseudogenes). These hits are tolerated, since
        % '(KRT17P1' begins with '(KRT17'.
        % 
        % However, it is possible that this may potentially introduce
        % artifacts. For example, when designing probes for CD4, you may
        % potentially have a hit in, say, CD44. In that case, since '(CD44'
        % begins with '(CD4', any such hit would be accepted. To check for
        % this potential artifact, the validation script
        % Targeted_ExSeq_ProbeGeneration_ValidateHomologyRegions.m can be
        % run to re-BLAST all of the identified homology regions for manual
        % verification that such errors did not occur. If unacceptable
        % sequences are identified, they can be removed for further probe
        % generation.
        %
        % Note that this behavior could potentially be modified if desired.
        % For example, when designing probes for ACTB (beta-actin), hits
        % within ACTA (alpha-actin) or ACTG (gamma-actin) genes could be
        % tolerated by modifying the code below.        
        
        matching_gene_names = strfind(blastout.Hits(ii).Name, ['(' gene_name]);
        matching_accession_numbers = strfind(blastout.Hits(ii).Name, accession_number);
                
        if ~isempty(matching_gene_names) || ~isempty(matching_accession_numbers)
            % Case: either the gene name or accession number of the gene of
            % interest is found in the hit. The query homology sequence is
            % not excluded.
                        
        else
            % Case: neither the gene name nor accession number of the gene
            % of interest are found in the hit, i.e. this is a hit in a
            % different gene. A more careful search is performed to
            % determine if the query sequence is excluded or not.
        
            for jj = 1:max(size(blastout.Hits(ii).HSPs))
                
                % Look for Plus/Plus hits, where the RNA strand segment
                % matches another RNA strand segment.
                is_plus_plus_stranded = strfind(blastout.Hits(ii).HSPs(jj).Strand, 'Plus/Plus');
                
                % Test to see if the number of matching bases of a hit is
                % greater or equal to than min_homology_length_to_exclude.
                has_matching_bases_hit = blastout.Hits(ii).HSPs(jj).Identities.Match >= min_homology_length_to_exclude;
                
                % Test for well-centered criteria: if the hit is
                % "well-centered", i.e. covers at least
                % ligation_junction_overlap_cutoff bases on both sides of
                % the ligation junction.
                is_well_centered = (fivePrime_end_position-blastout.Hits(ii).HSPs(jj).QueryIndices(1)) >= (ligation_junction_overlap_cutoff-1) &&  (blastout.Hits(ii).HSPs(jj).QueryIndices(2) - threePrime_start_position) >= ((ligation_junction_overlap_cutoff-1));
                               
                % If a hit is Plus/Plus, AND the hit has enough exactly
                % matching bases, AND the hit is well-centered, exclude the
                % query sequence. If all three conditions are not met, then
                % the homology region is acceptable.
                if is_plus_plus_stranded && has_matching_bases_hit && is_well_centered
                    is_excluded_for_homology = true;
                    
                    % optional lines that can be uncommented for
                    % visualization/debugging                             
                    %
                    % alignment = blastout.Hits(ii).HSPs(jj).Alignment
                    % blastout.Hits(ii).Name
                end
                
            end
        end
    end

    delete('temp.fa');    

end

