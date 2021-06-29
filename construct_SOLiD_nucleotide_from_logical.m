function [output_barcode_sequences] = construct_SOLiD_nucleotide_from_logical(SOLiD_logical_barcodes, barcode_specification, params)
    % construct_SOLiD_nucleotide_from_logical takes three arguments: (1) a 
    % cell array of SOlID logical barcodes (SOLiD_logical_barcodes) to
    % convert to nucleotide sequences; (2) a barcode specification
    % (barcode_specification) that defines the structure of the SOLiD
    % barcode, formatted as a table; (3) params, containing key parameters
    % for the function.
    %
    % The key parameter is solid_rca_binding_site
    % (params.RCA_PRIMER_BINDING_SITE).
    % 
    % The function returns a cell array containing the nucleotide sequences
    % that correspond to the logical barcodes.
    % 
    % The technical details of how this function works are explained
    % inline below.

    % set parameter    
    solid_rca_binding_site = params.RCA_PRIMER_BINDING_SITE;       
    
    % In the first step, we unpack the barcode specification (which is
    % described in terms of the order of ligations on specific primers) and
    % restructure it in terms of nucleotide positions.
             
    % Get sequencing primer names and ligation cycle names
    seq_primers = barcode_specification.Properties.RowNames;
    ligation_cycles = barcode_specification.Properties.VariableNames;
    
    % Get logical length of barcode from number of non-NaN values
    barcode_logical_length = sum(~isnan(barcode_specification{:,:}(:)));
           
    % Define barcode_indices, which is a cell array containing an unrolled
    % barcode specification. The rows correspond to position in the logical
    % barcode. The first column is the primer, the second column is the
    % ligation number, and the third column is the index of the first
    % position of the interrogated dibase.
    barcode_indices = cell(barcode_logical_length, 3); 
    
    % Unpack barcode_specification into barcode_indices. To do this, we
    % loop over barcode_specification, first by sequencing primer, then by
    % ligation number, and insert the appropriate values into the
    % appropriate location in barcode_indices.
    
    % Loop over sequencing primers
    for ii = 1:size(barcode_specification, 1)
        
        curr_seq_primer = seq_primers{ii};
        
        % Get index of number in curr_seq_primer (i.e. index of the 'X' in
        % 'SeqN-X'), and set curr_seq_primer_offset = X. 
        seq_primer_char = regexp(curr_seq_primer, '[0-9]');
        
        if ~isempty(seq_primer_char)
            curr_seq_primer_offset = str2num(curr_seq_primer(seq_primer_char));
        else
            curr_seq_primer_offset = 0; % for case of 'SeqN'
        end

        % Loop over ligations for the current sequencing primer
        for jj = 1:size(barcode_specification, 2)
            curr_barcode_position = barcode_specification{ii, jj};            
           
            % Record values in barcode_indices if it the barcode position
            % is utilized (i.e not a silent ligation).
            if ~isnan(curr_barcode_position)
                barcode_indices{curr_barcode_position, 1} = seq_primers{ii};
                barcode_indices{curr_barcode_position, 2} = ligation_cycles{jj};
                barcode_indices{curr_barcode_position, 3} = (jj-1)*5-curr_seq_primer_offset+1;
            end           
        end        
    end
    
    % Identify the starting position of the last ligation, which determines
    % overall barcode length
    max_starting_pos = max([barcode_indices{:,3}]);
    barcode_nucleotide_length = max_starting_pos + 1;    
    
    % Create cell array that will contain the final nucleotide sequences
    % for the logical barcodes passed to this function.
    output_barcode_sequences = cell(size(SOLiD_logical_barcodes, 1), 1);
    
    % Create a placeholder for the current barcode sequence to be computed.
    % We work in the space of the original padlock probe. We repeat 'N's
    % for the length of the barcode, and append the first four bases of the
    % RCA primer binding site, which is immedately 3' to the SOLiD barcode
    % on the padlock probe. We then compute the reverse complement, since
    % this is the template sequence in the amplicon that is decoded.
    padlock_sequence = [repmat('N', [1, barcode_nucleotide_length]), solid_rca_binding_site(1:4)];
    template_sequence = seqrcomplement(padlock_sequence); %0 position is at +4 in this array
       
    % Because SOLiD barcodes can potentially be degenerate if the coverage
    % is not complete (multiple nucleotide sequences corresponding to the
    % same colorspace barcode), we want to identify the barcode coverage to
    % ensure that it is unique. Non-unique barcodes are not a problem in
    % general, but they are not handled by this code. We compute the
    % barcode coverage by looping over the barcode indices.
    barcode_coverage = zeros(barcode_nucleotide_length+1, 1);
    
    all_starting_positions = [barcode_indices{:,3}]';
    
    for ii = 1:size(all_starting_positions, 1)
        barcode_coverage(all_starting_positions(ii)+1:all_starting_positions(ii)+2) = barcode_coverage(all_starting_positions(ii)+1:all_starting_positions(ii)+2) + 1;
    end
    barcode_coverage = barcode_coverage(2:end);
    
    if size(find(barcode_coverage == 0), 1) > 0
        error('Barcode not fully determined.');
    end
    
    % Now that we have established that the barcode coverage is complete,
    % we work through the barcode from close to the seq primer to far from
    % the seq primer, using the starting positions to determine which
    % logical barcode value to use to determine the current dibase.
    
    % Get the dibase start positions, and then store the index ordering of
    % lowest to highest start positions in min_start_ordering. This
    % ordering allows for the completion of the barcode from known
    % sequences (i.e. those that have an overlap with the seq primer)
    % outwards towards the unknown sequences.
    [min_values, min_start_ordering] = mink([barcode_indices{:,3}]', size(barcode_indices(:,3), 1));
    
    % Loop over each of the logical barcodes to translate
    for ii = 1:size(SOLiD_logical_barcodes,1)
        
        curr_SOLiD_bc = SOLiD_logical_barcodes{ii};
        
        % Strip leading quote from the barcode string
        if curr_SOLiD_bc(1) == ''''
            curr_SOLiD_bc = curr_SOLiD_bc(2:end);
        end
        
        curr_template = template_sequence;
        
        %disp(['Current barcode: ', curr_SOLiD_bc]);
        
        % Loop over the barcode positions from seq primer outwards to solve
        % the barcode sequence.        
        for jj = 1:size(min_start_ordering,1)
            
            % Get the current order position, then determine the logical
            % barcode value at the current position.
            curr_idx = min_start_ordering(jj);
            curr_logical_barcode = str2num(curr_SOLiD_bc(curr_idx));
            curr_nt_position = 4+barcode_indices{curr_idx,3};
            
            % The display commands can be uncommented to visualize the
            % barcode solution process.
            
            %disp(['Current logical barcode value: ', num2str(curr_logical_barcode)]);
            %disp(['Template nucleotide sequence: ', curr_template]);            
            %disp(['Position: ' num2str(4+barcode_indices{curr_idx,3})]);
            %disp(['Dibase: ' num2str(curr_template(4+barcode_indices{curr_idx,3}:4+barcode_indices{curr_idx,3}+1))]);            
            
            % Use the mapping given by SOLiD_logical_to_dibase to solve the
            % the identity of the current dibase sequence under evaluation.
            % Since solving the nucleotide sequence proceeds from the
            % sequencing primer out, each dibase under consideration has
            % one known nucleotide sequence and one unknown nucleotide
            % sequence. In each iteration, we identify the unknown
            % nucleotide sequence. In the next iteration, we shift the
            % dibsae window by one nucleotide, so the recently-identified
            % nucleotide now becomes the known nucleotide for the next
            % iteration.
            %
            % SOLiD_logical_to_dibase takes the current logical barcode
            % value at the current position, and the current dibase
            % sequence, and returns the final sequence for that dibase.
            
            curr_template(curr_nt_position:curr_nt_position+1) = SOLiD_logical_to_dibase(curr_logical_barcode, curr_template(curr_nt_position:curr_nt_position+1), params);
                        
            %disp(['Updated sequence: ', curr_template]);                                                                  
        end
        
        % Strip the sequencing primer and compute reverse complement of the
        % sequence to translate back from templatespace to the sequence to 
        % insert on the padlock probe.
        output_barcode_sequences{ii} = seqrcomplement(curr_template(5:end));        
    end
end