function [Name, WLevel, WPID] = DWNodeName(Level, Index)
% function [Name, WLevel, WPID] = DWNodeName(Level, Index)
%
% Return the symbolic name of a node in a diffusion wavelet tree at the
% given level and index.
%
% In:
%    Level = level
%    Index = node index
% Out:
%    Name  = a string indicating the symbolic name of the node; e.g. 'W_3'
%
% Dependences:
%    none

% check to see if index is valid
if Index <= 2^ceil(Level/2)

   WLevel = [];
   WPID = [];

   if Index == 1
      Name = sprintf('V_%d', Level);
   elseif Index == 2
      Name = sprintf('W_%d', Level);
   else

      % we have a wavelet packet node
      % find the level of its parent W
      WLevel = Level;
      WPID = '';

      while Index~=2
         WLevel = WLevel-1;
         if mod(Index,2)==0
            WPID = [WPID '1'];
         else
            WPID = [WPID '0'];
         end
         Index=ceil(Index/2);
      end
      Name = sprintf('W_%d (%s)', WLevel, WPID);
   end
end


