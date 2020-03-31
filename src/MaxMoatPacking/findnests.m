%findsubsets
%should also sort, prior to finding

function [mAdj, M] = findnests(w, M)
wAdj = zeros(length(w));
counter = 0;
modi = zeros(length(w));
[ht, wdth] = size(M);
mAdj = zeros(ht, wdth);
nestInList = zeros(wdth, wdth+1);
for i=1:wdth
    for j=i+1:wdth
      %  Mnew = M
        potential=0;
        iSubsetofj = 0;
        jSubsetofi = 0;
        subs = 1;
   
        for k=1:ht
            if (M(k, i) + M(k,j) == 0)
                potential = 1;
            end
        end
        
        if (potential==0) 
            continue;
        end
        for k=1:ht
            
            if (M(k, i) > M(k,j))
                iSubsetofj = 1;
            end
            if (M(k,i) < M(k,j))
                jSubsetofi = 1;
            end
        end
        
     
            
        
        if (iSubsetofj == 1)
            for k=1:ht
                 if (M(k, i) < M(k,j))
                     subs=0;
                     break;
                 end
            end
            wAdj(i) = wAdj(i)+wAdj(j);
        end
        if (jSubsetofi == 1)
               for k=1:ht
                 if (M(k, i) > M(k,j))
                     subs=0;
                     break;
                 end
               end
            if (subs == 1)
                if (modi(i) == 1)
                    break;
                end
                modi(i) = 1;
                subsetRows =  M(:,11)+M(:,14)==0;
                if (iSubsetofj == 1)
                    mAdj(:,i) = mAdj(:,i) + subsetRows * w(j);
                    currNestedIn = nestInList(j,1);
                    nestInList(j, currNestedIn + 2) = i;
                    nestInList(j, 1) = currNestedIn + 1;
                end
                if (jSubsetofi == 1)
                    mAdj(:,j) = mAdj(:,j) + subsetRows * w(i);
                    currNestedIn = nestInList(i,1);
                    nestInList(i, currNestedIn + 2) = j;
                     nestInList(i, 1) = currNestedIn + 1;
                end
                
                counter = counter + 1;
           %     Mnew(:, i) = M(:, j);
            %    Mnew(:, j) = M(:, i);
            end
        end
    end
 %  M = Mnew;

end
 %disp(counter)
  rem = [];
 size(nestInList);
 for v = 1:wdth+1
     if (sum(nestInList(:, v)) == 0)
         rem = [rem, v];
     end
 end
 wMod = zeros(ht,1);
 nestInList(:, rem) = [];
 disp(nestInList)
 for i=1:wdth
           if (nestInList(i, 1) == 0)
          continue;
      end
    for j=1:max(1,nestInList(i, 1))
        curr = nestInList(i, j);
        if (curr == 0) 
            break;
        end
        wMod(curr, 1) = wMod(curr, 1) + w(i);
    end
 end
 del = zeros(wdth);
 wdthF = wdth;
  for i=1:wdth
      if (nestInList(i, 1) == 0)
          continue;
      end
    for j=1:nestInList(i, 1)
        curr = nestInList(i, j);
        if (i > curr)
          
            M(:, wdthF+1) = M(:, j);
            wdthF = wdthF + 1;
            del(j) = del(j) + 1;
        end
    end
  end
 Mholder = M;

  for i=1:wdth
      if (del(i) > 0)
          M(i) = [];
          del(i) = del(i) - 1;
      end
  end
 w = wMod
 disp(ht)
 