function Prop2time(Q,no_region)
for i = 1:no_region
    x = sum(abs(Q(i,:)));
    if x >= 1
        flag = 1
    end
end
end