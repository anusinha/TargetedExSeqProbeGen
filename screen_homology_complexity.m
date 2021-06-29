function [excluded_for_complexity] = screen_homology_complexity(query_sequence)
    % screen_homology_complexity takes a query sequence (query_sequence)
    % and returns a boolean variable excluded_for_complexity that
    % represents if the query sequence passed the minimum sequence
    % complexity requirements.
    % 
    % The sequence complexity requirements are that there are no repeats of
    % length 5 or longer of the same nucleotide, and that all four
    % nucleotides are represented in the sequence.

    excluded_for_complexity = false;
   
    if seqwordcount(query_sequence,'AAAAA') > 0 | seqwordcount(query_sequence,'CCCCC') > 0 | seqwordcount(query_sequence,'GGGGG') > 0 | seqwordcount(query_sequence,'TTTTT') > 0
        excluded_for_complexity = true;
    end

    number_of_different_nucleotides = logical(seqwordcount(query_sequence,'T')) + logical(seqwordcount(query_sequence,'A')) + ...
        logical(seqwordcount(query_sequence,'C')) + logical(seqwordcount(query_sequence,'G'));
    if number_of_different_nucleotides < 4
        excluded_for_complexity = true;
    end
end

