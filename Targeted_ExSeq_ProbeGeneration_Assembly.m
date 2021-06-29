% This script is the second of two to run to generate barcoded padlock
% probes for a multiplexed Targeted ExSeq experiment. This script assumes
% that the the genelist has been constructed, containing the barcode
% assignment for each transcript, and that the first script has been run,
% which constructs the homology regions for padlock probes.
%
% This script converts the barcode sequences (designed in a logical barcode
% space) to nucleotide sequences, and assembles the final padlock probes.
% The script also outputs the padlock probes sequences in formats that are
% suitable for uploading to IDT for ordering in plates.
%
% For further details, see the probe generation protocol.

close all; clear;

%% START OF USER INPUT -- SECTION 1 %%

% Required parameters to set:
% 1) species
% 2) metadata
% 3) date_of_generation
% 4) genelist_spreadsheet
% 5) homology_regions_directory
% 6) output_directory
% 7) all_probes_output_filename
% 8) probes_output_by_sheets_filename
% 9) plate_map_for_pooling_filename

% Optional parameters to change (all parameters saved in params struct):
% 1) RCA_PRIMER_BINDING_SITE = 'TCTCGGGAACGCTGAAGACGGC'
% 2) SPACER_SEQUENCE = 'AAA'
% 3) MAX_NUMBER_OF_PROBES = 16
% 4) PLATE_NUM_ROWS = 8
% 5) PLATE_NUM_COLS = 12
% 6) AUTOGENERATE_BARCODES = true

% Computed parameters
% 1) SOLID_SEQN_PRIMER (from RCA_PRIMER_BINDING_SITE)
% 2) PLATE_SIZE (from PLATE_NUM_ROWS x PLATE_NUM_COLS)

%%%% REQUIRED PARAMETERS TO SET %%%%

% Summary of probe generation: species, type of experiment (metadata), and
% date of generation are recorded for plate record information
species = ''; % i.e. Human, Mouse, etc
metadata = ''; % i.e. Brain, BreastCancer, etc 
date_of_generation = ''; % i.e. 20210101
    
% Paths to input files

% Filename for genelist Excel file. This Excel file has four columns:
% gene symbol, RefSeq accession number, minimum Tm threshold, and barcode.
% Each row corresponds to a gene with an assigned barcode. See protocol for
% futher details
genelist_spreadsheet = '';

% Name of directory containing the possible binding (homology) regions
% (the output of the first script) for each gene in the genelist.
homology_regions_directory = '';    

% Paths to output files

% Name of the directory in which the Excel files with the probes will be
% generated.
output_directory = '';

% Name of the output Excel file containing all of the probes (as one long
% spreadsheet).
all_probes_output_filename = '';

% Name of the output Excel file containing all of the probes organized into
% plates for upload to IDT.
probes_output_by_sheets_filename = '';

% Name of the output Excel file containing maps of all of the plates, to
% help guide pooling probes.
plate_map_for_pooling_filename = '';

%%%% OPTIONAL PARAMETERS TO MODIFY %%%%

% PROBE GENERATION PARAMETERS

% This sequence is on the backbone of the padlock probe, and is the
% sequence that the rolling circle amplification primer will hybridize to.
% The default sequence (TCTCGGGAACGCTGAAGACGGC) is a modified version of
% the original universal primer from FISSEQ to increase its Tm, improving
% its compatibility with Illumina sequencing. This sequence should not have
% any spaces in it. The SOLiD SeqN primer is a subset of this sequence, and
% is computed below, saved in SOLID_SEQN_PRIMER.
params.RCA_PRIMER_BINDING_SITE = 'TCTCGGGAACGCTGAAGACGGC';
params.SOLID_SEQN_PRIMER = params.RCA_PRIMER_BINDING_SITE(1:18);

% This spacer sequence is added between the homology arms of the padlock
% probes and the SOLiD/Illumina barcodes (added 5' to the SOLiD barcode, 3'
% to the Illumina barcode). The default sequence (AAA) is selected to
% set the total length of the padlock probe to be 74 nt.
params.SPACER_SEQUENCE = 'AAA';

% The maximum number of padlock probes to generate per gene.
params.MAX_NUMBER_OF_PROBES = 16;

% PLATE OUTPUT PARAMETERS

% Define the number of rows and number of columns of the plate.
params.PLATE_NUM_ROWS = 8;
params.PLATE_NUM_COLS = 12;

% A boolean variable representing whether or not to automatically generate
% the SOLiD and Illumina barcodes (as described by the protocol).
params.AUTOGENERATE_BARCODES = true;

% Computed parameters
params.PLATE_SIZE = params.PLATE_NUM_ROWS * params.PLATE_NUM_COLS;

%%%% END OF USER INPUT -- SECTION 1 %%%%
%% START OF USER INPUT -- SECTION 2 %%

%%%% ADVANCED PARAMETERS TO MODIFY %%%%
%
% These parameters govern the nature of the mapping between logical barcode
% values and nucleotide bases for the SOLiD and Illumina sequencing
% chemistries. These parameters generally do not need to be changed if the
% 4-color Illumina chemistry or SOLiD sequencing chemistries are being used
% as-is.

% Optional advanced parameters to change:
% A) Parameters relevant to Illumina sequencing
%     1) ILLUMINA_MAP_LOGICAL
%     2) ILLUMINA_MAP_TEMPLATE
% B) Parameters relevant to SOLiD sequencing
%     1) SOLID_NUCLEOTIDE_KEYS
%     2) SOLID_DIBASE_MATRIX
%     
% Computed parameters
% 1) ILLUMINA_MAP_PADLOCK (from ILLUMINA_MAP_TEMPLATE)
% 2) ILLUMINA_MAP (stored in params struct) (from ILLUMINA_MAP_LOGICAL and ILLUMINA_MAP_PADLOCK)
% 3) SOLID_DIBASE_ENCODING_SPECIFICATION (stored in params struct, from SOLID_NUCLEOTIDE_KEYS and SOLID_DIBASE_MATRIX)
% 4) SOLID_LOGICAL_KEYS (from SOLID_DIBASE_ENCODING_SPECIFICATION)
% 5) SOLID_LOGICAL2DIBASE_VALS (from SOLID_LOGICAL_KEYS and SOLID_DIBASE_ENCODING_SPECIFICATION)
% 6) SOLID_LOGICAL2DIBASE (stored in params struct, from SOLID_LOGICAL_KEYS and SOLID_LOGICAL2DIBASE_VALS)

% Illumina Sequencing Barcoding Specification
% 
% We implement a direct mapping between logical barcodes and nucleotides
% for the Illumina chemistry. The mapping between logical barcodes and
% nucleotides is defined by the variables ILLUMINA_MAP_LOGICAL and
% ILLUMINA_MAP_TEMPLATE, where the template refers to the sequencing
% template, which is the rolling circle amplification product. We compute
% the padlock nucleotides in ILLUMINA_MAP_PADLOCK based on the template
% nucleotides, and construct the ILLUMINA_MAP which is used by the
% function construct_Illumina_nucleotide_from_logical.
%
% The standard mapping maps the logical bases {'0', '1', '2', '3'} to the
% padlock bases {'A', 'C', 'G', 'T'} (with the template bases being {'T',
% 'G', 'C', 'A'}.
%
% See protocol for further details on the mapping between the logical
% barcode space and the Illumina chemistry.   

ILLUMINA_MAP_LOGICAL = {'0', '1', '2', '3'};
ILLUMINA_MAP_TEMPLATE = {'T', 'G', 'C', 'A'};

ILLUMINA_MAP_PADLOCK = cellstr(cellfun(@seqcomplement, ILLUMINA_MAP_TEMPLATE'))';

params.ILLUMINA_MAP = containers.Map(ILLUMINA_MAP_LOGICAL, ILLUMINA_MAP_PADLOCK);

% SOLiD Sequencing Barcoding Specification
%
% For the SOLiD sequencing chemistry, we implement a direct mapping between
% dibase nucleotides and logical barcodes. The reverse mapping between
% logical barcodes and dibases is degenerate, as there are multiple dibases
% that correspond to a single logical barcode value. This mapping is used
% by other functions that translate a logical barcode into a nucleotide
% sequence.
%
% The mapping between dibase nucleotides and logical barcodes is shown in
% the table below. The first template base (of the dibase) determines the
% row, and the second template base (of the dibase) determines the column.
% For example, for the template dibase 5' AC 3', the first base of the
% dibase is A, and the second is C, corresponding to row 1, column 2 in the
% table, which is logical 1.
%
%           2nd base
%           A  C  G  T
%         +-----------
%       A | 0  1  2  3
% 1st   C | 1  0  3  2
% base  G | 2  3  0  1
%       T | 3  2  1  0
%
% Numbers are the logical barcode value and correspond to physical colors
% (in conventional SOLiD sequencing) as:
%    0 - 488 nm excitation
%    1 - 561 nm
%    2 - 594 nm
%    3 - 647 nm
% i.e. 0 - 3 corresponds to lowest to highest wavelength.
%
% This scheme is implemented with the following variables. The nucleotides
% are defined by SOLID_NUCLEOTIDE_KEYS, and the logical barcode matrix is
% defined numerically in SOLID_DIBASE_MATRIX.
% 
% The table is assembled in SOLID_DIBASE_ENCODING_SPECIFICATION (stored in
% the params struct so that it can be passed to other functions). 
% 
% Logical keys are identified in SOLID_LOGICAL_KEYS, and the reverse
% mapping from logical keys to dibases is computed in the map
% SOLID_LOGICAL2DIBASE (stored in the params struct).

% Define the nucleotides for the row/columns of the SOLiD matrix, and 
% matrix values.
SOLID_NUCLEOTIDE_KEYS = {'A'; 'C'; 'G'; 'T'};
SOLID_DIBASE_MATRIX = [[0, 1, 2, 3]; [1, 0, 3, 2]; [2, 3, 0, 1]; [3, 2, 1, 0]];

% Construct SOLID_DIBASE_ENCODING_SPECIFICATION 
params.SOLID_DIBASE_ENCODING_SPECIFICATION = array2table(SOLID_DIBASE_MATRIX);
params.SOLID_DIBASE_ENCODING_SPECIFICATION.Properties.RowNames = SOLID_NUCLEOTIDE_KEYS;
params.SOLID_DIBASE_ENCODING_SPECIFICATION.Properties.VariableNames = SOLID_NUCLEOTIDE_KEYS;

% Identify logical keys from table representation
SOLID_LOGICAL_KEYS = num2cell(unique(params.SOLID_DIBASE_ENCODING_SPECIFICATION{:,:}(:)));

% Find dibase nucleotides corresponding to logical barcode values
SOLID_LOGICAL2DIBASE_VALS = cell(size(SOLID_LOGICAL_KEYS, 1), 1);  
for ii = 1:size(SOLID_LOGICAL_KEYS, 1)
    [row_idx, col_idx] = ind2sub([size(params.SOLID_DIBASE_ENCODING_SPECIFICATION, 1), size(params.SOLID_DIBASE_ENCODING_SPECIFICATION, 2)], find(table2array(params.SOLID_DIBASE_ENCODING_SPECIFICATION) == SOLID_LOGICAL_KEYS{ii}));
    row_names = params.SOLID_DIBASE_ENCODING_SPECIFICATION.Properties.RowNames(row_idx);
    col_names = params.SOLID_DIBASE_ENCODING_SPECIFICATION.Properties.VariableNames(col_idx);
    SOLID_LOGICAL2DIBASE_VALS{ii} = cellfun( @(x,y) [x, y], row_names, col_names', 'UniformOutput', false)';        
end

% Create map between logical keys and dibase cell arrays
params.SOLID_LOGICAL2DIBASE = containers.Map(SOLID_LOGICAL_KEYS, SOLID_LOGICAL2DIBASE_VALS);

%%%% END OF USER INPUT -- SECTION 2 %%%%
%% START OF USER INPUT -- SECTION 3 %%

%%%% SOLID SEQUENCING-SPECIFIC BARCODE SPECIFICATION PARAMETERS %%%%

% Because SOLiD sequencing uses dibase interrogation and sequences bases in
% the barcode in a non-sequential order, the nucleotide sequence of the
% barcode depends on how the order of ligations on different primers
% corresponds to different positions in the logical barcode. We define the
% correspondence between physical bases of sequencing and logical bases
% below.
% 
% A barcode library of length 7 and minimum Hamming Distance 3 enables
% error detection and one-bit error correction for ~340 targets. We
% implement this barcoding approach with SOLiD in two potential ways,
% although the end user can change the barcode specification to suit 
% their need.
%
% Barcode Readout Approach 1: efficient length-4 subset readout (as
% highlighted in protocol)
%
% In the first approach, the first four logical barcode positions
% correspond to the first two ligations on SeqN-1, followed by the first
% two ligations on SeqN. The last three positions correspond to the second
% ligations on SeqN-2, N-3, and N-4, respectively. This is shown
% schematically in the figure below. This utility of this barcode
% specification is that a smaller (~50 barcodes) length 4, Hamming
% Distance 2 barcode library can be embedded in the first four bases of the
% longer library. This enables the user to start with a smaller gene set of
% ~50 genes that can be read out in four bases of sequencing for a pilot
% experiment. The user can then add additional targets (up to ~340),
% that can be read out in seven rounds of sequencing with imaging
% (requiring 10 ligations, with three of the ligations being silent
% (non-imaged) because of SOLiD sequencing.
%
% The advantage of this barcoding scheme is that it allows the subset of
% length-4 barcodes to be read out in four rounds of sequencing, without
% any silent ligations. 
%
% This implementation is shown below. The rows indicate the seq primers,
% and the columns indicate bases along the template (5' -> 3', with bases
% 1-7 corresponding to the seven nucelotide bases of the barcode.
% Interrogated dibases are shown with pairs of numbers in adjacent bases,
% with the number corresponding to the position in the logical barcode.
% Silent ligations are indicated with parentheses around the barcode number
% they immediately preceed, i.e. to access barcode position 5 (SeqN-2,
% Ligation 2), you first need to pass through the silent position (5)
% (SeqN-2, Ligation 1). x's denote positions that are not identifiable with
% that seq primer.
% 
%    Table for Approach 1: efficient length-4 subset readout
%
%               -3  -2  -1   0   1   2   3   4   5   6   7
%            +--------------------------------------------
%        N   |                   3   3   x   x   x   4   4  
%        N-1 |               1   1   x   x   x   2   2   x
%        N-2 |          (5) (5)  x   x   x   5   5   x   x
%        N-3 |      (6) (6)  x   x   x   6   6   x   x   x
%        N-4 |  (7) (7)  x   x   x   7   7   x   x   x
% 
% Barcode Readout Approach 2: consecutive embedding of bases (alternative
% approach mentioned in protocol)
%
% This approach is conceptually simpler: consecutive dibases correspond to
% the next position of the SOLiD barcode, i.e. dibase 0:1 is the logical
% position 1, dibase 1:2 is the logical position 2, and so on up to dibase
% 6:7 being logical position 7. This approach is conceptually simpler, and
% is included because we used it previously. 
% 
% Note that the experimental protocol for reading out the full length
% barcode is the same, as the first two ligations on each of the seq
% primers need to be performed. However, the second ligations on SeqN-1 and
% SeqN correspond to logical bases 6 and 7, respectively. For reading out
% the embedded length 4 library, it requires two silent ligations to be
% performed (on SeqN-4 and SeqN-3, respectively), so this scheme is less
% efficient. 
%
% While length-4 barcodes could potentially be embedded here, we suggest
% using approach 1 for ease of designing an embedded barcode set, and for
% the sake of simplicity and consistency across Illumina and SOLiD
% barcodes. The barcode design is shown diagramatically below (same
% (notation as above).
%
%    Table for Approach 2: consecutive embedding of bases
%    
%               -3  -2  -1   0   1   2   3   4   5   6   7
%            +--------------------------------------------
%        N   |                   2   2   x   x   x   7   7  
%        N-1 |               1   1   x   x   x   6   6   x
%        N-2 |          (5) (5)  x   x   x   5   5   x   x
%        N-3 |      (4) (4)  x   x   x   4   4   x   x   x
%        N-4 |  (3) (3)  x   x   x   3   3   x   x   x
%

% Optional advanced parameters to change to change SOLiD barcode
% specification:
% 1) MAX_NUM_LIGATIONS (in params)
% 2) The actual barcoding scheme (see note below)
% 3) T_selected_barcode_specification
% 4) barcode_descriptor


% The maximum number of successive ligations on a single sequencing primer
% needed to implement the barcoding scheme.
params.MAX_NUM_LIGATIONS = 2;

% Define primer names (rows) and ligation labels (columns) for barcoding
% tables.
LIGATION_PRIMERS = {'SeqN', 'SeqN-1', 'SeqN-2', 'SeqN-3', 'SeqN-4'};
LIGATION_POSITIONS = strcat('Ligation_', cellstr(num2str([1:params.MAX_NUM_LIGATIONS]'))');

% Create a row vector of NaNs used to define a table of NaNs (before
% specifying values).
variableTypes = repmat({'doublenan'}, 1, params.MAX_NUM_LIGATIONS);

% Approach 1: efficient length-4 subset readout (as described in the
% protocol)
T_barcode_specification_embedded_len4_HD2 = table('Size', [5, params.MAX_NUM_LIGATIONS], 'VariableTypes', variableTypes, 'RowNames', LIGATION_PRIMERS, 'VariableNames', LIGATION_POSITIONS);
T_barcode_specification_embedded_len4_HD2{'SeqN-1', 'Ligation_1'} = 1;
T_barcode_specification_embedded_len4_HD2{'SeqN-1', 'Ligation_2'} = 2;
T_barcode_specification_embedded_len4_HD2{'SeqN',   'Ligation_1'} = 3;
T_barcode_specification_embedded_len4_HD2{'SeqN',   'Ligation_2'} = 4;
T_barcode_specification_embedded_len4_HD2{'SeqN-2', 'Ligation_2'} = 5;
T_barcode_specification_embedded_len4_HD2{'SeqN-3', 'Ligation_2'} = 6;
T_barcode_specification_embedded_len4_HD2{'SeqN-4', 'Ligation_2'} = 7;

% Approach 2: consecutive embedding of bases (alternative approach)
T_barcode_specification_consecutive_bases = table('Size', [5, params.MAX_NUM_LIGATIONS], 'VariableTypes', variableTypes, 'RowNames', LIGATION_PRIMERS, 'VariableNames', LIGATION_POSITIONS);
T_barcode_specification_consecutive_bases{'SeqN-1', 'Ligation_1'} = 1;
T_barcode_specification_consecutive_bases{'SeqN',   'Ligation_1'} = 2;
T_barcode_specification_consecutive_bases{'SeqN-4', 'Ligation_2'} = 3;
T_barcode_specification_consecutive_bases{'SeqN-3', 'Ligation_2'} = 4;
T_barcode_specification_consecutive_bases{'SeqN-2', 'Ligation_2'} = 5;
T_barcode_specification_consecutive_bases{'SeqN-1', 'Ligation_2'} = 6;
T_barcode_specification_consecutive_bases{'SeqN',   'Ligation_2'} = 7;

% Users can specify their own barcoding schemes here.


% Select the barcoding scheme used in the remainder of the code.
T_selected_barcode_specification = T_barcode_specification_embedded_len4_HD2

% Add additional metadata to the full probe name about the SOLiD barcoding
% specification used.
solid_barcode_descriptor = 'embedded_len4';

%%%% END OF USER INPUT -- SECTION 3 %%%%
%%
% Save current warning state and temporarily suppress fastawrite warnings
warnState = warning;
warning('off','Bioinfo:fastawrite:AppendToFile');

% Check that all mandatory parameters are specified
if isempty(species) || isempty(metadata) || isempty(date_of_generation) ...
        || isempty(genelist_spreadsheet) || isempty(output_directory) ...
        || isempty(all_probes_output_filename) || isempty(probes_output_by_sheets_filename) ...
        || isempty(plate_map_for_pooling_filename)
    error('Did not specify mandatory parameter.');
end

% Load and unpack genelist and barcode library
T_genelist = readtable(genelist_spreadsheet, 'ReadVariableNames', false, 'ReadRowNames', false);
T_genelist.Properties.VariableNames = {'Symbol', 'Accession', 'MinTm', 'Barcode'};
all_gene_names = T_genelist.Symbol;
all_accession_numbers = T_genelist.Accession;

% Unpack barcodes
all_barcodes = T_genelist.Barcode;
all_barcodes = cellfun(@(x) x(2:end), all_barcodes, 'UniformOutput', false); % strip leading quote

num_genes = size(T_genelist, 1);

% Create output directory
mkdir(output_directory);

% Create cell that will contain gene names and a cell containing probe
% names and probe sequences for each gene
C_generated_probes = cell(num_genes, 2);

% Loop over each gene and generate probes for each gene
for ii = 1:num_genes
    tic
    % Get curent gene name, accession number, barcode
    gene_name = all_gene_names{ii};
    accession = strtok(all_accession_numbers{ii}, '.'); % strtok to strip characters from the additional info in the accession beyond the . point
    barcode = all_barcodes{ii};       
        
    % Option to automatically generate logical SOLiD and Illumina barcodes
    % given the barcoding scheme and the barcode of interest
    if params.AUTOGENERATE_BARCODES == true
        % Generate SOLiD Barcode for gene of interest
        SOLiDBarcode = construct_SOLiD_nucleotide_from_logical({barcode}, T_selected_barcode_specification, params);
        SOLiDBarcode = SOLiDBarcode{:};

        % Generate Illumina Barcode for gene of interest
        IlluminaBarcode = construct_Illumina_nucleotide_from_logical({barcode}, params);
        IlluminaBarcode = IlluminaBarcode{:};
    else
        % Option for user to do something else with non-autogenerated
        % barcodes (i.e. using custom code to generate barcodes).
        error('No manual generation for barcode sequences specified.');
    end

    % Read in homology regions (generated in previous script) for current
    % gene of interest
    try
        all_homology_regions = fastaread(['./' homology_regions_directory '/' accession '.fa']);
    catch        
        warning(['no read probes from accession: ' accession]);
        all_homology_regions = [];
    end

    % To generate probes, loop over all the possible binding sites for each
    % gene and generate all of the potential probes, up to the limit of
    % MAX_NUMBER_OF_PROBES. Probes are generated from homology regions
    % starting at 5' end of transcript, proceeding towards 3' end.
    
    % Create a cell with dimensions of MAX_NUMBER_OF_PROBES x 3 to store the
    % full probe name, short probe name, and probe sequence for each
    % generated probe (with empty string cells for probes that have not
    % been generated).
    curr_gene_probe_output = cellfun(@(x) '', cell(params.MAX_NUMBER_OF_PROBES,3), 'UniformOutput', false);

    % Loop over the homology regions, up to the limit of
    % MAX_NUMBER_OF_PROBES, and generate probe names and sequences, saving
    % the results in curr_gene_probe_output.
    for jj = 1:min([length(all_homology_regions) params.MAX_NUMBER_OF_PROBES])       
        
        % Get homology sequence, get length, and identify 5' and 3' arms
        homology_region = all_homology_regions(jj).Sequence;
        
        homology_length = size(homology_region, 2);
        threePrime_end = floor(homology_length/2);
        fivePrime_start = threePrime_end + 1;
        
        % Assemble full probe name:
        % species_metadata_date_of_generation_geneSymbol_accession_design_solidBarcodeDescriptor_barcode_#######_probe_#
        % Assemble short probe name:
        % species_metadata_symbol_bc_#######_probe_#
        full_probe_name = [species, '_', metadata, '_', date_of_generation, '_', gene_name, '_', accession, '_design_', solid_barcode_descriptor '_barcode_', barcode, '_probe_' num2str(jj)];
        short_probe_name = [species, '_', metadata, '_', gene_name, '_bc_', barcode, '_probe_' num2str(jj)];
        
        % Standard barcode autogeneration
        if params.AUTOGENERATE_BARCODES == true
            % probe assembly: 5' phosphorylation, 5' homology arm, spacer,
            % SOLiD barcode, RCA primer binding site, Illumina barcode,
            % spacer, 3' homology arm
            probe = ['/5Phos/', homology_region(fivePrime_start:homology_length), params.SPACER_SEQUENCE, SOLiDBarcode, params.RCA_PRIMER_BINDING_SITE, IlluminaBarcode, params.SPACER_SEQUENCE, homology_region(1:threePrime_end)];
        else
            error('No custom probe architecture for manual barcode generation specified');
        end
        
        % Save probe output
        curr_gene_probe_output{jj, 1} = full_probe_name;
        curr_gene_probe_output{jj, 2} = short_probe_name;
        curr_gene_probe_output{jj, 3} = probe;

    end
    % Save gene name and cell array containing probe names and probe
    % sequences in C_generated_probes.
    C_generated_probes{ii, 1} = gene_name;
    C_generated_probes{ii, 2} = curr_gene_probe_output;
end

% Prior warning state is restored
warning(warnState); %Reset warning state to previous settings

%% Export probes to Excel spreadsheets for output

% Exported spreadsheets:
% 1) All probes output (all_probes_output_filename) -- lists all probe
% sequences in a single spreadsheet.
% 2) Probes output by sheets (probes_output_by_sheets_filename) -- outputs
% probes into individual 96-well plate sheets, formatted for upload to IDT
% for synthesis in plate format
% 3) Plate map for pooling (plate_map_for_pooling_filename) -- outputs
% maps of plates to help guide the user in pooling plates of probes
% together.

% (1) All probes spreadsheet: this spreadsheet lists all probes in a single
% Excel sheet, one probe per row. The columns are: (1) gene symbol; (2)
% full probe name; (3) short probe name; (4) probe sequence. If the number
% of probes generated for a transcript is below MAX_NUMBER_OF_PROBES, empty
% rows are included in the spreadsheet.

% Pull data into C_all_probes cell and write cell to spreadsheet file
C_all_probes = cell(num_genes * params.MAX_NUMBER_OF_PROBES, 4);
for ii = 1:num_genes
    for jj = 1:params.MAX_NUMBER_OF_PROBES
        
        C_all_probes{(ii-1)*params.MAX_NUMBER_OF_PROBES + jj, 2} = C_generated_probes{ii, 2}{jj, 1};
        C_all_probes{(ii-1)*params.MAX_NUMBER_OF_PROBES + jj, 3} = C_generated_probes{ii, 2}{jj, 2};
        C_all_probes{(ii-1)*params.MAX_NUMBER_OF_PROBES + jj, 4} = C_generated_probes{ii, 2}{jj, 3};
        
        if ~isempty(C_generated_probes{ii, 2}{jj, 1})
            C_all_probes{(ii-1)*params.MAX_NUMBER_OF_PROBES + jj, 1} = C_generated_probes{ii, 1};
        end
    end
end

writecell(C_all_probes, ['./' output_directory '/' all_probes_output_filename])

% (2) Probes output by sheets and (3) Plate map for pooling spreadsheets
% are generated together. A requirement is that the number of rows on the
% the plates (8) is a divisor of MAX_NUMBER_OF_PROBES.

% Pad C_all_probes with empty rows so that the total number of rows is
% divisible by PLATE_SIZE. This allows the cell array to be reshaped.

total_number_of_plates = ceil(size(C_all_probes, 1)/params.PLATE_SIZE);
rows_to_pad = total_number_of_plates * params.PLATE_SIZE - size(C_all_probes, 1);

C_all_probes_padded = [C_all_probes; repmat({''}, rows_to_pad, 4)];

% Extract specific features, and reshape as 3D matrices with dimensions:
% PLATE_NUM_ROWS x PLATE_NUM_COLS x total_number_of_plates.
C_all_probes_genes_reshaped = reshape(C_all_probes_padded(:,1), params.PLATE_NUM_ROWS, params.PLATE_NUM_COLS, total_number_of_plates);
C_all_probes_full_names_reshaped = reshape(C_all_probes_padded(:,2), params.PLATE_NUM_ROWS, params.PLATE_NUM_COLS, total_number_of_plates);
C_all_probes_short_names_reshaped = reshape(C_all_probes_padded(:,3), params.PLATE_NUM_ROWS, params.PLATE_NUM_COLS, total_number_of_plates);
C_all_probes_sequences_reshaped = reshape(C_all_probes_padded(:,4), params.PLATE_NUM_ROWS, params.PLATE_NUM_COLS, total_number_of_plates);

% Generate annotation matrix
plate_map_row_labels = cellstr(char('A' + [0:params.PLATE_NUM_ROWS-1]'))';
plate_map_col_labels = strtrim(cellstr(num2str([1:params.PLATE_NUM_COLS]'))');

% single_plate_labels contains labels A1, A2,...; B1, B2,...;...H11, H12 in
% plate dimensions.
single_plate_labels = cell(params.PLATE_NUM_ROWS, params.PLATE_NUM_COLS);
for ii = 1:params.PLATE_NUM_ROWS
    for jj = 1:params.PLATE_NUM_COLS
        well_row = char('A' + (ii-1));
        well_col = num2str(jj);        
        single_plate_labels{ii,jj} = [well_row well_col];        
    end
end
single_plate_labels_unrolled = single_plate_labels(:);

% Write probe tables: loop over the total number of plates, generate the
% output for each plate, then write out the current sheet, before moving
% onto the next plate.
for ii = 1:total_number_of_plates    
    % Pull out features for current plate being processed
    curr_plate_gene_names = C_all_probes_genes_reshaped(:,:,ii);    
    curr_plate_full_names = C_all_probes_full_names_reshaped(:,:,ii);
    curr_plate_full_names = curr_plate_full_names(:);
    curr_plate_probe_sequences = C_all_probes_sequences_reshaped(:,:,ii);
    curr_plate_probe_sequences = curr_plate_probe_sequences(:);
      
    % Assemble current plate containing probe sequences (in IDT's upload
    % format)
    T_curr_plate = table(single_plate_labels_unrolled, curr_plate_full_names, curr_plate_probe_sequences, 'VariableNames', {'Well Position', 'Name', 'Sequence'});
    
    % Assemble plate maps containing gene names in occupied wells
    T_plate_annotations = array2table(curr_plate_gene_names, 'RowNames', plate_map_row_labels, 'VariableNames', plate_map_col_labels);
    
    % Write out to specific sheets in spreadsheet
    writetable(T_curr_plate, ['./', output_directory, './', probes_output_by_sheets_filename], 'WriteVariableNames', true, 'WriteRowNames', false, 'Sheet', ['Plate_', num2str(ii), '_', species, '_', metadata]);
    writetable(T_plate_annotations, ['./', output_directory, './', plate_map_for_pooling_filename], 'WriteVariableNames', true, 'WriteRowNames', true, 'Sheet', ['Plate_', num2str(ii), '_', species, '_', metadata]);    
end

disp('Output Probe Spreadsheets Written');