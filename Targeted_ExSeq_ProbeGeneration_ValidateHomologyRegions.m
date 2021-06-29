% This script is an optional script to validate the homology regions
% identified by the Targeted_ExSeq_ProbeGeneration_Homology.m script.
% This script assumes the existence of a genelist, relevant RefSeq
% database, and a directory containing the homology region output of the
% previous script.
% 
% For further details, see the probe generation protocol.

close all; clear;

%% START OF USER INPUT %%

% Required parameters to set:
% 1) genelist_spreadsheet
% 2) RefSeq_database
% 3) homology_regions_output

% Constants set (optional parameters to modify, all within the params
% struct; values should exactly match the parameters previously used in
% Targeted_ExSeq_ProbeGeneration_Homology):
% 1) HOMOLOGY_LENGTH = 32
% 2) E_VALUE = 0.05
% 3) BLAST_HOMOLOGY_MIN_LENGTH_CUTOFF = 13
% 4) BLAST_HOMOLOGY_LIGATION_JUNCTION_MIN_OVERLAP_CUTOFF = 4

%%%% REQUIRED PARAMETERS TO SET %%%%

% Filename for genelist Excel file. This Excel file has four columns:
% gene symbol, RefSeq accession number, minimum Tm threshold, and barcode.
% Each row corresponds to a gene with an assigned barcode. See protocol for
% futher details
genelist_spreadsheet = '';

% Filename for the reference RefSeq database for BLAST. Resulting files
% (.fna, .nhr, .nin, .nsq files) should be saved in same folder as this
% script. All files should have same filename before file extension. See
% protocol for further details.
refseq_database = '';

% Input directory of homology region sequences (output of
% Targeted_ExSeq_ProbeGeneration_Homology.m) where the possible homology
% regions for each gene in the genelist are saved. Note that the sequences
% in this directory are the reverse complements of the original RNA
% sequences.
homology_regions_input = '';

% Output filename for spreadsheet sumarizing homology region validation
output_filename = '';

%%%% CONSTANTS (OPTIONAL PARAMETERS TO MODIFY) %%%%

% Default length of the homology region is 32, corresponding to 16 bases on
% each arm of the padlock probe.
params.HOMOLOGY_LENGTH = 32;

% The maximum acceptable E value for BLAST queries to be considered for
% further evaluation. For a particular hit, a higher E value indicates a
% less-specific match; lower E values correspond to more specific matches.
% This parameter sets the threshold E value for which hits are retained for
% further evaluation. See protocol for further note on E values.
params.E_VALUE = 0.05;

% In the homology search, the minimum number of matching bases of a BLAST
% hit to proceed with considering query sequence exclusion. BLAST hits
% with fewer matching bases are accepted.
params.BLAST_HOMOLOGY_MIN_LENGTH_CUTOFF = 13;

% In the homology search, the minimum overlap of the BLAST hit on both
% sides of the ligation junction to exclude the query. 
params.BLAST_HOMOLOGY_LIGATION_JUNCTION_MIN_OVERLAP_CUTOFF = 4;

% Computed parameters
params.FIVEPRIME_END = floor(params.HOMOLOGY_LENGTH/2);
params.THREEPRIME_START = params.FIVEPRIME_END + 1;

%%%% END OF USER INPUT %%%%
%%
% Save current warning state
warnState = warning;
warning('off', 'Bioinfo:fastawrite:AppendToFile');
warning('off', 'MATLAB:DELETE:FileNotFound');

% Check that all mandatory parameters are specified
if isempty(genelist_spreadsheet) || isempty(refseq_database) || isempty(homology_regions_input)
    error('Did not specify mandatory parameter.');
end

% Load genelist table
T_genelist = readtable(genelist_spreadsheet, 'ReadVariableNames', false, 'ReadRowNames', false);
T_genelist.Properties.VariableNames = {'gene', 'accession', 'minTm', 'barcode'};

% Set up summary output cell array. The output for each gene is vertically
% concatenated to this cell array. Column 1 = gene name; column 2 =
% homology sequence (reverse complement of RNA transcript); column 3 = cell
% array of gene names of all hits.

summary_output = cell(0,3);

% Loop over all genes, and BLAST all homology regions
for genes = 1:size(T_genelist,1)
tic    
    gene_name = T_genelist.gene{genes};
    accession = T_genelist.accession{genes};
    accession = strtok(accession, '.');
    
    % gene_specific_output cell: this cell will contain the BLAST results
    % for the homology regions for the current gene of interest.
    % We will add a row for each homology region. Column 1 = gene name;
    % column 2 = homology sequence (reverse complement of RNA transcript);
    % column 3 = cell array of gene names of all hits.
    gene_specific_output = cell(0, 3);    
    
    % Read homology regions
    try
        all_homology_regions = fastaread(['./' homology_regions_input '/' accession '.fa']);        
    catch
        warning('no read probes from accession');
        all_homology_regions = [];
    end       
    
    disp(['Running validation for gene ' num2str(genes) ' of ' num2str(size(T_genelist,1)) ': ' gene_name]);    
    
    % Loop over homology regions for current gene of interest
    for sequence = 1:size(all_homology_regions, 1)
        
        % Reverse complement of the saved sequence so that we are BLASTing
        % the actual RNA sequence
        query_sequence = seqrcomplement(all_homology_regions(sequence).Sequence);
        gene_hits = {};
        
        delete('temp.fa');            
        fastawrite('temp.fa', query_sequence); 
    
        blastlocal_parameters = ['-S 1 -e ', num2str(params.E_VALUE)];
    
        % Run BLAST query and save output in blastout. This assumes that the
        % BLAST has been added to the PATH variable. 
        blastout = blastlocal('InputQuery','temp.fa','Program','blastn','DataBase', refseq_database, 'blastargs', blastlocal_parameters); 
        
        for hits = 1:size(blastout.Hits, 2)
            currentHit = blastout.Hits(hits);
            
            for alignment = 1:length(currentHit.HSPs)
                currentAlignment = blastout.Hits(hits).HSPs(alignment);
                
                % Look for Plus/Plus hits, where the RNA strand segment
                % matches another RNA strand segment.
                is_plus_plus_stranded = strfind(currentAlignment.Strand, 'Plus/Plus');
                                
                % Test to see if the number of matching bases of a hit is
                % greater or equal to than min_homology_length_to_exclude.
                has_matching_bases_hit = currentAlignment.Identities.Match >= params.BLAST_HOMOLOGY_MIN_LENGTH_CUTOFF;
                
                % Test for well-centered criteria: if the hit is
                % "well-centered", i.e. covers at least
                % ligation_junction_overlap_cutoff bases on both sides of
                % the ligation junction.
                is_well_centered = (params.FIVEPRIME_END-currentAlignment.QueryIndices(1)) >= (params.BLAST_HOMOLOGY_LIGATION_JUNCTION_MIN_OVERLAP_CUTOFF-1) &&  (currentAlignment.QueryIndices(2) - params.THREEPRIME_START) >= ((params.BLAST_HOMOLOGY_LIGATION_JUNCTION_MIN_OVERLAP_CUTOFF-1));
                               
                if is_plus_plus_stranded && has_matching_bases_hit && is_well_centered
                    tokens = split(blastout.Hits(hits).Name, {'(', ')'});
                    hit_gene_name = tokens{2};
                    gene_hits{1,end+1} = hit_gene_name;
                end
            end
        end
        
        % Deduplication of gene hits in the same gene (i.e. if a sequence
        % matches multiple variants)
        gene_hits = unique(gene_hits);
        % format output for current homology region.
        currSequenceOutput = {gene_name, all_homology_regions(sequence).Sequence, gene_hits};
        gene_specific_output = [gene_specific_output; currSequenceOutput];
        delete('temp.fa');
    end
    
    summary_output = [summary_output; gene_specific_output];   
toc    
end

% We now unroll the summary_output cell array so we can write it to an
% Excel spreadsheet.
summary_output_unrolled = cell(0,0);

for ii = 1:size(summary_output, 1)
    summary_output_unrolled(ii, 1) = summary_output(ii, 1);
    summary_output_unrolled(ii, 2) = summary_output(ii, 2);    
    numHits = size(summary_output{ii,3}, 2);
    for jj = 1:numHits
        if ~isempty(summary_output{ii,3})        
            summary_output_unrolled{ii,2+jj} = summary_output{ii,3}{jj};
        end
    end    
end

writecell(summary_output_unrolled, output_filename);

warning(warnState); % Reset warning state to previous settings
