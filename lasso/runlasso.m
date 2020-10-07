#! /usr/bin/octave-cli -q

## Prefix for the ACP files. Leave empty if you do not want the ACP files written. 
prefix="lasso";

## List of 1-norm constraints to use
tlist = linspace(0,0.3,21);
#tlist = [0.1 0.15 0.3 0.5 1 2 3 5 10];

## Use the maximum coefficients, if available?
usemaxcoef=0;

#### Do NOT touch past here ####
if (!exist("octavedump.m","file"))
  error("This script requires an octavedump.m file to work.")
endif
source("octavedump.m");

## scale the ACP terms with the coef0
for i = 1:columns(x)
  x(:,i) /= coef0(i);
endfor

## start the loop
nacp = 0;
printf("#it    t           norm-1       norm-2      norm-inf      wrms    nterm iter\n");
for it = 1:length(tlist)
  t = tlist(it);

  ## scale the columns using the maximum coefficient information
  if (exist("maxcoef","var") && !isempty(maxcoef) && t > 0 && exist("usemaxcoef","var") && usemaxcoef)
    for i = 1:columns(x)
      xtilde(:,i) = x(:,i) * maxcoef(i) / mean(maxcoef);
    endfor
    [w,iteration] = LassoActiveSet(xtilde,y,t);
    for i = 1:columns(x)
      w(i) = w(i) * maxcoef(i) / mean(maxcoef);
    endfor
    clear xtilde
  else
    [w,iteration] = LassoActiveSet(x,y,t,'optTol',1e-9,'zeroThreshold',1e-9);
  endif

  worig = w;

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
    for i = 1:length(w)
      fprintf(fid,"%.15e\n",w(i));
    endfor
    fclose(fid);
  endif

  ## Calculate the wrms and print the line
  wrms = sqrt(sum((y - x * worig).^2));
  printf("%2.2d %.10f %.10f %.10f %.10f %.10f %d %d\n",it,t,...
         sum(abs(w)),sqrt(sum(w.^2)),max(abs(w)),...
         wrms,sum(sum(nterms)),iteration);
endfor

