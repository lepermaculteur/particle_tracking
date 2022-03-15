function insert_disk(pos,r,col)
  
  if strcmp(col,'blue')
    plot(pos(:,1),pos(:,2),'o','MarkerFaceColor','b', ... 
      'MarkerEdgeColor','b','MarkerSize',r)  
  endif

  if strcmp(col,'red')
    plot(pos(:,1),pos(:,2),'o','MarkerFaceColor','r', ... 
      'MarkerEdgeColor','r','MarkerSize',r)    
  endif

endfunction

##  [XX, YY] = meshgrid (1:bord, bord:-1:1);
##
##  for k=1:length(pos(:,1))
##    mask = ((XX-pos(k,1)).^2 + (YY-pos(k,2)).^2) <= r^2;
##
##    if (strcmp(col,'red'))
##      mask = cat (3, mask, false (size (mask)), false (size (mask)));
##    elseif (strcmp(col,'blue'))
##      mask = cat (3, false (size (mask)), false (size (mask)), mask);
##    endif
##
##    im(mask)=255;
##  endfor  
##  if strcmp(col,'blue')   
##    
##    im(max(round(pos(1,2))-1,1):min(round(pos(1,2))+1,bord),max(round(pos(1,1))-3,1),3)=256;
##    im(max(round(pos(1,2))-2,1):min(round(pos(1,2))+2,bord),max(round(pos(1,1))-2,1),3)=256;
##    im(max(round(pos(1,2))-3,1):min(round(pos(1,2))+3,bord),max(round(pos(1,1))-1,1),3)=256;
##    im(max(round(pos(1,2))-3,1):min(round(pos(1,2))+3,bord),round(pos(1,1)),3)=256;
##    im(max(round(pos(1,2))-3,1):min(round(pos(1,2))+3,bord),min(round(pos(1,1))+1,bord),3)=256;
##    im(max(round(pos(1,2))-2,1):min(round(pos(1,2))+2,bord),min(round(pos(1,1))+2,bord),3)=256;
##    im(max(round(pos(1,2))-1,1):min(round(pos(1,2))+1,bord),min(round(pos(1,1))+3,bord),3)=256;
##    
##    
##    im(max(round(pos(1,2))-1,1):min(round(pos(1,2))+1,bord),max(round(pos(1,1))-3,1),1:2)=0;
##    im(max(round(pos(1,2))-2,1):min(round(pos(1,2))+2,bord),max(round(pos(1,1))-2,1),1:2)=0;
##    im(max(round(pos(1,2))-3,1):min(round(pos(1,2))+3,bord),max(round(pos(1,1))-1,1),1:2)=0;
##    im(max(round(pos(1,2))-3,1):min(round(pos(1,2))+3,bord),round(pos(1,1)),1:2)=0;
##    im(max(round(pos(1,2))-3,1):min(round(pos(1,2))+3,bord),min(round(pos(1,1))+1,bord),1:2)=0;
##    im(max(round(pos(1,2))-2,1):min(round(pos(1,2))+2,bord),min(round(pos(1,1))+2,bord),1:2)=0;
##    im(max(round(pos(1,2))-1,1):min(round(pos(1,2))+1,bord),min(round(pos(1,1))+3,bord),1:2)=0;
##    
##  elseif strcmp(col,'red')
##    
##    im(max(round(pos(1,2))-1,1):min(round(pos(1,2))+1,bord),max(round(pos(1,1))-3,1),1)=256;
##    im(max(round(pos(1,2))-2,1):min(round(pos(1,2))+2,bord),max(round(pos(1,1))-2,1),1)=256;
##    im(max(round(pos(1,2))-3,1):min(round(pos(1,2))+3,bord),max(round(pos(1,1))-1,1),1)=256;
##    im(max(round(pos(1,2))-3,1):min(round(pos(1,2))+3,bord),round(pos(1,1)),1)=256;
##    im(max(round(pos(1,2))-3,1):min(round(pos(1,2))+3,bord),min(round(pos(1,1))+1,bord),1)=256;
##    im(max(round(pos(1,2))-2,1):min(round(pos(1,2))+2,bord),min(round(pos(1,1))+2,bord),1)=256;
##    im(max(round(pos(1,2))-1,1):min(round(pos(1,2))+1,bord),min(round(pos(1,1))+3,bord),1)=256;
##    
##    
##    im(max(round(pos(1,2))-1,1):min(round(pos(1,2))+1,bord),max(round(pos(1,1))-3,1),2:3)=0;
##    im(max(round(pos(1,2))-2,1):min(round(pos(1,2))+2,bord),max(round(pos(1,1))-2,1),2:3)=0;
##    im(max(round(pos(1,2))-3,1):min(round(pos(1,2))+3,bord),max(round(pos(1,1))-1,1),2:3)=0;
##    im(max(round(pos(1,2))-3,1):min(round(pos(1,2))+3,bord),round(pos(1,1)),2:3)=0;
##    im(max(round(pos(1,2))-3,1):min(round(pos(1,2))+3,bord),min(round(pos(1,1))+1,bord),2:3)=0;
##    im(max(round(pos(1,2))-2,1):min(round(pos(1,2))+2,bord),min(round(pos(1,1))+2,bord),2:3)=0;
##    im(max(round(pos(1,2))-1,1):min(round(pos(1,2))+1,bord),min(round(pos(1,1))+3,bord),2:3)=0;
##    
##  endif
