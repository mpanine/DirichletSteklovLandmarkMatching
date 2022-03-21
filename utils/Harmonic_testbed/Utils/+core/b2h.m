function [ hex_str ] = b2h( bin_vec )
    %BINARYVECTORTOHEX Converts a binary vector to a hex string.
    
    if mod(numel(bin_vec),4) ~= 0
        throw(MException('BINARYVECTORTOHEXERROR:badInput',...
            'The input of binaryVectorToHex has to be a multiple of 4.'));
    end
    
%     fprintf('B2H');
    hex_str = b2h_core(bin_vec);
%     hex_str = '';
%     for i=1:4:numel(bin_vec)
%         b = strrep(num2str(bin_vec(i:i+3)),' ','');
%         h = '';
%         switch b
%             case '0000'
%                 h = '0';
%             case '0001'
%                 h = '1';
%             case '0010'
%                 h = '2';
%             case '0011'
%                 h = '3';
%             case '0100'
%                 h = '4';
%             case '0101'
%                 h = '5';
%             case '0110'
%                 h = '6';
%             case '0111'
%                 h = '7';
%             case '1000'
%                 h = '8';
%             case '1001'
%                 h = '9';
%             case '1010'
%                 h = 'A';
%             case '1011'
%                 h = 'B';
%             case '1100'
%                 h = 'C';
%             case '1101'
%                 h = 'D';
%             case '1110'
%                 h = 'E';
%             case '1111'
%                 h = 'F';
% %             case '1010'
% %                 h = 'a';
% %             case '1011'
% %                 h = 'b';
% %             case '1100'
% %                 h = 'c';
% %             case '1101'
% %                 h = 'd';
% %             case '1110'
% %                 h = 'e';
% %             case '1111'
% %                 h = 'f';
%         end
%         hex_str = [hex_str,h];
%     end
end

