% minMax2 determines the minimum and maximum of all data in the two input% arrays. Minmax(1) is the minimum, minmax(2) is the maximum.% minMax2 is for instance useful for producing equal axes on x-y graphsfunction minmax = minMax2(data1,data2)minmax(1) = min([data1(:); data2(:)]);minmax(2) = max([data1(:); data2(:)]);