function [output_nucleotides] = construct_Illumina_nucleotide_from_logical(Illumina_logical_barcodes, params)
    % construct_Illumina_nucleotide_from_logical takes a cell array of
    % logical barcodes, and the parameters struct, containing the mapping
    % between logical barcodes and nucleotides, and returns a cell array of
    % nucleotide barcodes, output_nucleotides, that represents the logical
    % barcodes in the space of Illumina sequencing.
    %
    % The output nucleotide sequence is directly placed into the padlock
    % probe, immediately 3' to the RCA primer binding site.
    
    output_nucleotides = Illumina_logical_barcodes;
    
    for key = params.ILLUMINA_MAP.keys
        keyChar = char(key);
        output_nucleotides = replace(output_nucleotides, key, params.ILLUMINA_MAP(keyChar));
    end
end