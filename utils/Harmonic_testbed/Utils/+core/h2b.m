function [ bin_vec ] = h2b( hex_str, n_bits )
    %HEXTOBINARYVECTOR Converts hex string to binary vector.
%     fprintf('H2B');

    bin_vec = h2b_core(hex_str, n_bits);
    
%     if numel(hex_str) < n_bits/4
%         missing_num = n_bits/4 - numel(hex_str);
%         hex_str = [repmat('0', 1, missing_num), hex_str];
%     elseif numel(hex_str) > n_bits/4
%         throw(MException('HEXTOBINARYVECTORERROR:badInput',...
%             'The input of hexToBinaryVector has to be a hex string with n_bits/4 characters.'));
%     end
%     
%     bin_vec = zeros(1, n_bits);
%     
%     for i=1:numel(hex_str)
%         h = hex_str(i);
%         b = '';
%         switch h
%             case {'0'}
%                 b = '0000';
%             case {'1'}
%                 b = '0001';
%             case {'2'}
%                 b = '0010';
%             case {'3'}
%                 b = '0011';
%             case {'4'}
%                 b = '0100';
%             case {'5'}
%                 b = '0101';
%             case {'6'}
%                 b = '0110';
%             case {'7'}
%                 b = '0111';
%             case {'8'}
%                 b = '1000';
%             case {'9'}
%                 b = '1001';
%             case {'A', 'a'}
%                 b = '1010';
%             case {'B', 'b'}
%                 b = '1011';
%             case {'C', 'c'}
%                 b = '1100';
%             case {'D', 'd'}
%                 b = '1101';
%             case {'E', 'e'}
%                 b = '1110';
%             case {'F', 'f'}
%                 b = '1111';
%         end
%         
%         bin_vec(1, 4*(i-1)+1:4*(i-1)+4) = arrayfun(@(x)str2double(x),b);
%     end
    
end

