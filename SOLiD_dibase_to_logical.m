function [outputLogical] = SOLiD_dibase_to_logical(inputDibase, solid_dibase_encoding_specification)
    % Given a nucleotide dibase and a table representing the SOLiD
    % sequencing dibase encoding specification, the SOLiD_dibase_to_logical
    % function returns the correponding logical colorspace value for SOLiD
    % sequencing. The mapping between template bases and colors is defined
    % by the following matrix (for standard SOLiD sequencing).
    %
    %     A  C  G  T
    %   +-----------
    % A | 0  1  2  3
    % C | 1  0  3  2
    % G | 2  3  0  1
    % T | 3  2  1  0
    %
    % Rows indicate the first base of the dibase (closer to the sequencing 
    % primer); second base of the dibase is in the columns.    
    %
    % The mapping can also be recast in terms of the logical barcode value
    % for each template dibase as follows:
    % 0: AA, CC, GG, TT
    % 1: AC, CA, GT, TG
    % 2: AG, CT, GA, TC
    % 3: AT, CG, GC, TA
    %
    % Numbers are the logical barcode value and correspond to physical
    % colors (in conventional SOLiD sequencing) as:
    %    0 - 488 nm excitation
    %    1 - 561 nm
    %    2 - 594 nm
    %    3 - 647 nm
    % i.e. 0 - 3 corresponds to lowest to highest wavelength.
    
    outputLogical = solid_dibase_encoding_specification(inputDibase(1), inputDibase(2));

end