function [excluded_for_physical_properties] = screen_homology_physical(query_sequence, minTmThreshold, params)
    % screen_homology_physical takes a query sequence (query_sequence), and
    % the parameters struct, and returns a boolean variable 
    % excluded_for_physical_properties that represents if the query
    % sequence passed the physical property sequence requirements.
    %
    % The relevant parameters used are:
    % (1) fivePrime_end and threePrime_start, which define the endpoints of
    % the independent ends of the padlock probe;
    % (2) gc_range, a two element array, where gc_range(1) is the minumum
    % GC content, and gc_range(2) is the maximum GC content;
    % (3) minTmThreshold and maxTmSeparation, which define the minimum Tm
    % theshold and minumum Tm separation between the two arms.    
    % 
    % The physical property requirements are in three parts: (1) GC range
    % for each individual arm; (2) Tm for each individual arm; and (3)
    % hairpin/dimer secondary structure for the entire sequence. See
    % further details below.    

    fivePrime_end = params.FIVEPRIME_END;
    threePrime_start = params.THREEPRIME_START;
    gc_range = params.GC_RANGE;
    maxTmSeparation = params.MAX_TM_SEPARATION;
        
    excluded_for_physical_properties = false;

    % Generate the two individual arms of the padlock probe binding site
    fivePrimeArm = query_sequence(1:fivePrime_end);
    threePrimeArm = query_sequence(threePrime_start:end);
    
    % Compute physical properties of the individual 5' and 3' arms,
    % assuming a salt concentration of 0.3M, and a padlock probe
    % concentration of 100 nM (per probe).
    fivePrimeArmSeqProperties = oligoprop(fivePrimeArm, 'Salt', 0.3, 'Primerconc', 100e-9);
    threePrimeArmSeqProperties = oligoprop(threePrimeArm, 'Salt', 0.3, 'Primerconc', 100e-9);

    % Screening GC content criteria: test that the GC content of each
    % individual arm is within the range specified in gc_range   
    if fivePrimeArmSeqProperties.GC < gc_range(1) || gc_range(2) < fivePrimeArmSeqProperties.GC
        excluded_for_physical_properties = true;
    end
    if threePrimeArmSeqProperties.GC < gc_range(1) || gc_range(2) < threePrimeArmSeqProperties.GC
        excluded_for_physical_properties = true;
    end
    
    % Screening melting temperature criteria:
    %   (1) test for melting temperatures greater than the Tm threshold
    % (minTmThreshold) for each individual padlock probe arm;
    %   (2) test that the melting temperatures for each individual padlock
    %   probe arm are within the minimum Tm separation (maxTmSeparation)
    %
    % The Tm is computed in the oligoprop function above. Note that this Tm
    % refers to the DNA-DNA Tm, and the final binding interaction for the
    % padlock probe to the transcript is DNA-RNA, which will be slightly
    % higher. We use the oligoprop Tm values 2 through 6, which use
    % salt-adjusted and nearest-neighbor models to compute the melting
    % temperature.
    
    fivePrimeArmTm = mean(fivePrimeArmSeqProperties.Tm(2:6));
    threePrimeArmTm = mean(threePrimeArmSeqProperties.Tm(2:6));
    
    if fivePrimeArmTm > minTmThreshold && threePrimeArmTm > minTmThreshold && abs(fivePrimeArmTm - threePrimeArmTm) <= maxTmSeparation
        
    else
        excluded_for_physical_properties = true;
    end

    % Screening secondary stucture criteria: no hairpins or dimers. Minimum
    % length of hairpin base is 7, with hairpin loop length of 6, and
    % minimum dimer length is 10.
    fullSeqProperties = oligoprop(query_sequence, 'HPBase', 7, 'HPLoop', 6, 'Dimerlength', 10);  
    if ~isempty(fullSeqProperties.Hairpins) || ~isempty(fullSeqProperties.Dimers)
        excluded_for_physical_properties = true;
    end

end

