distance=[60,65,70,75,80,85,90,95,100];
cnt=5;
s1=zeros(9,1);
s=zeros(cnt,1);
for i=1
    for j=1:cnt
        s(j)=main_dist(distance(i)); 
    end
    s1(i)=mean(s);
end
