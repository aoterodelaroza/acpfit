#! /usr/bin/octave-cli -q

## Prefix for the ACP files. Leave empty if you do not want the ACP files written. 
prefix="lasso";

## List of 1-norm constraints to use
tlist = linspace(0,10,21);
#tlist = [0.1 0.15 0.3 0.5 1 2 3 5 10];
#tlist = [0 4];
#tlist = [100];

## Use the maximum coefficients, if available?
usemaxcoef=0;

#### Do NOT touch past here ####

## the version of this lasso script
lasso_version = "1.2";

## Read the dump binary file if availble. If not, read the text file and generate the binary file.
filedump = "octavedump.m";
finfo = dir(filedump);
if (isempty(finfo))
  error("This script requires an octavedump.m file to work.")
endif
filebin = sprintf("%s.bin",filedump);
fbininfo = dir(filebin);
if (!isempty(fbininfo) && fbininfo.statinfo.mtime >= finfo.statinfo.mtime)
  load("-binary",filebin);
else
  source(filedump);
  save("-binary",filebin,"atoms","lmax","lname","explist","nrows","ncols","coef0","x","y","maxcoef",...
      "yaddnames","yadd");
endif

## place the additional energy terms as columns
if (exist("yadd","var"))
  x = [x, yadd];
endif

## start the loop
nacp = 0;
printf("| #id |   t      |     norm-1   |    norm-2   |   norm-inf   |   wrms  |  nterm | iter |\n");
for it = 1:length(tlist)
  t = tlist(it);

  ## scale the columns using the maximum coefficient information
  if (exist("maxcoef","var") && !isempty(maxcoef) && t > 0 && exist("usemaxcoef","var") && usemaxcoef)
    factor = ones(columns(x),1);
    xtilde = zeros(size(x));
    maxcoef = maxcoef(:);

    while (1)
      for i = 1:columns(x)
        xtilde(:,i) = x(:,i) * factor(i);
      endfor
      [w,iteration] = LassoActiveSet(xtilde,y,t,'optTol',1e-9,'zeroThreshold',1e-9);
      for i = 1:length(w)
        w(i) = w(i) * factor(i);
      endfor

      wrms = sqrt(sum((y - x * w).^2));
      idx = find(abs(w) > maxcoef);
      printf("norm-1 = %.4f, wrms = %.4f, maxcoef: ok = %d / violations = %d\n",sum(abs(w)),wrms,columns(x)-length(idx),length(idx));
      if (length(idx) == 0)
        break
      endif

      for i = 1:length(w)
        if (abs(w(i)) > maxcoef(i)) 
          ratio = min(maxcoef(i) / abs(w(i)),0.99);
          factor(i) = factor(i) * ratio^(1.1);
        endif
      endfor
    endwhile
    clear xtilde factor changed
  else
    [w,iteration] = LassoActiveSet(x,y,t,'optTol',1e-9,'zeroThreshold',1e-9);
  endif

  ## Calculate the wrms and print the line
  wrms = sqrt(sum((y - x * w).^2));

  ## unpack the ACP coefficients from the row vector
  nterms = zeros(length(atoms),length(lname));
  n = 0;
  for iat = 1:length(atoms)
    for il = 1:lmax(iat)
      for iexp = 1:length(explist)
        n++;
        coef = w(n);
        if (abs(coef) > 1e-20)
          nterms(iat,il)++;
          aexp(iat,il,nterms(iat,il)) = explist(iexp);
          acoef(iat,il,nterms(iat,il)) = coef;
        endif
      endfor
    endfor
  endfor

  ## Write the ACP in Gaussian format
  if (exist("prefix","var") && !isempty(prefix))
    nacp++;
    fid = fopen(sprintf("%s-%2.2d.acp",prefix,nacp),"w");

    ## Write the ACP header
    fprintf(fid,"! This ACP was generated with runlasso.m, version %s.\n",lasso_version);
    fprintf(fid,"! Atoms: ");
    for i = 1:length(atoms)
      fprintf(fid,"%2.2s ",atoms{i});
    endfor
    fprintf(fid,"\n");
    fprintf(fid,"! Lmax:  ");
    for i = 1:length(lmax)
      fprintf(fid,"%2.2s ",lname{lmax(i)});
    endfor
    fprintf(fid,"\n");
    fprintf(fid,"! Exponents: ");
    for i = 1:length(explist)
      fprintf(fid,"%.2f ",explist(i));
    endfor
    fprintf(fid,"\n");
    fprintf(fid,"! ACP terms in training set: %d\n",ncols);
    fprintf(fid,"! Data points in training set: %d\n",nrows);
    if (exist("yadd","var"))
      fprintf(fid,"! Additional energy contributions:\n");
      for i = 1:length(yaddnames)
        fprintf(fid,"!  %s: coeff = %.10f\n",yaddnames{i},w(ncols+i));
      endfor
    endif
    if (exist("maxcoef","var") && !isempty(maxcoef) && t > 0 && exist("usemaxcoef","var") && usemaxcoef)
      fprintf(fid,"! Maximum coefficients applied\n");
    endif
    fprintf(fid,"! norm-1 = %.4f\n",sum(abs(w)));
    fprintf(fid,"! norm-inf = %.4f\n",max(abs(w)));
    fprintf(fid,"! wrms = %.4f\n",wrms);

    ## Write the ACP itself
    for i = 1:length(atoms)
      fprintf(fid,"%s 0\n",atoms{i});
      fprintf(fid,"%s %d 0\n",atoms{i},lmax(i)-1);
      for j = 1:lmax(i)
        fprintf(fid,"%s\n",lname{j});
        fprintf(fid,"%d\n",nterms(i,j));
        for k = 1:nterms(i,j)
          fprintf(fid,"2 %.6f %.15f\n",aexp(i,j,k),acoef(i,j,k));
        endfor
      endfor
    endfor
    fclose(fid);

    fid = fopen(sprintf("%s-%2.2d.w",prefix,nacp),"w");
    for i = 1:ncols
      fprintf(fid,"%.15e\n",w(i));
    endfor
    fclose(fid);
  endif

  printf("| %2.2d |  %.4f |  %.4f | %.4f | %.4f | %.4f | %d | %d |\n",it,t,...
         sum(abs(w)),sqrt(sum(w.^2)),max(abs(w)),...
         wrms,sum(sum(nterms)),iteration);
endfor

