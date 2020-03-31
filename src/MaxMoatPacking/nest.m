% DE NESTING TEST CASE

%w = [1, 2, 3, 4, 0]';
%mvt = [1 1 1 1 1;
%      0 1 0 1 1;
 %     0 1 1 0 1;
 %      0 0 1 1 1;
 %      1 1 1 1 1;
 %      1 0 1 1 0;
 %      1 0 1 1 0];

 function [wo, mvto] = nest(w, mvt)
  
w = w';
cont = 1;
while (cont ==1)
  [ht, wdth] = size(mvt);
  if (wdth <= 1) 
      break;
  end
  change = 0;
  for i=1:wdth-1
      
      
      lencomp = 0;
      rencomp = 0;
      overlap = 0;
      for j=1:ht
            if ((mvt(j,i) == 0) && (mvt(j,i+1)==1))
                lencomp = 1;
            end
            if ((mvt(j,i) == 1) && (mvt(j,i+1)==0))
                rencomp = 1;
            end
             if ((mvt(j,i) == 0) && (mvt(j,i+1)==0))
                overlap = 1;
            end
      end
      if (lencomp + rencomp + overlap == 3)
         epsilon = min(w(i), w(i+1));
         
         removeMoat = i;
         if (w(i+1) < w(i))
             removeMoat = i+1;
         end
         
         w(i) = w(i) - epsilon;
         w(i+1) = w(i+1) - epsilon;

         S1 = ones(ht, 1) - (mvt(:,i) < mvt(:, i+1));
         S2 = ones(ht, 1) - (mvt(:,i) > mvt(:, i+1));
         wS1 = epsilon;
         wS2 = epsilon;
         
         S1new = 1;
           
         for k=1:wdth
          
             
             if (sum(S1-mvt(:,k) == 0) == ht)
                 w(k) = w(k) + wS1;
                 S1new = 0;
                 break;
             end
     
         end
         
         S2new = 1;
         for k=1:wdth
             

             if (sum(S2-mvt(:,k) == 0) == ht)
                 w(k) = w(k) + wS2;
                 S2new = 0;
                 break;
             end
     
         end
         
         if (S1new == 1)
             mvt = [mvt, S1];
             w = [w, wS1];
             change = 1;
             
         end
         if (S2new == 1)
             mvt = [mvt, S2];
             w = [w, wS2];
             change = 1;
          
         end
            w(removeMoat) = [];
            mvt(:, removeMoat) = [];
            [ht, wdth] = size(mvt);
            change = 1;
            break;
        end  
  end
  if (change == 0)
  cont = 0;
  end
 end
  wo = w';
  mvto=mvt;
  
  
  