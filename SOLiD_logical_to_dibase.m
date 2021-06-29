function [output_nucleotides] = SOLiD_logical_to_dibase(logicalBC, diBase_template, params)
    % Given a logical barcode value and an optional dibase template, this
    % function returns the dibase(s) that could generate the logical
    % barcode. This function requires the params struct; the value used is
    % the SOLID_LOGICAL2DIBASE Map.
    %
    % If one argument is provide (just the logical barcode value), all four
    % dibases corresponding to that logical barcode value are returned. 
    % If two arguments are provided, the logical barcode value and the
    % dibase template (of the form 'XN' or 'NX', where 'X' is one of
    % {'A', 'C', 'G', 'T'} that specifies a fixed/constrained base), the
    % one matching dibase corresponding to that logical value is returned.
    % If 'NN' is used as the template, four dibases are returned. 
    %
    % The mapping between template bases and colors is defined by the
    % SOLiD dibase encoding specification. The standard SOLiD logical
    % barcodes are shown in the following following matrix.
    %
    %     A  C  G  T
    %   +-----------
    % A | 0  1  2  3
    % C | 1  0  3  2
    % G | 2  3  0  1
    % T | 3  2  1  0
    %
    % Rows indicate the first base of the dibase (closer to the sequencing 
    % primer); second base of the dibase is in the colums.    
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
    
    switch nargin
        case 2
            output_nucleotides = params.SOLID_LOGICAL2DIBASE(logicalBC);
        case 3
            logical_dibases = params.SOLID_LOGICAL2DIBASE(logicalBC);
            
            if strcmp(diBase_template, 'NN')
                output_nucleotides = params.SOLID_LOGICAL2DIBASE(logicalBC);
            elseif diBase_template(1) == 'N'
                output_nucleotides = logical_dibases{diBase_template(2) == cellfun(@(x) x(2), logical_dibases')};
            elseif diBase_template(2) == 'N'
                output_nucleotides = logical_dibases{diBase_template(1) == cellfun(@(x) x(1), logical_dibases')};              
            end
    end
end

