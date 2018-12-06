tl=[];
for j=2:1500
    t=r_rec(:,j)'*r_rec(:,j-1);
    tl=[tl,t];
end