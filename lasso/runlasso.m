#! /usr/bin/octave-cli -q

## Prefix for the ACP files. Leave empty if you do not want the ACP files written. 
prefix="lasso";

## List of 1-norm constraints to use
tlist = linspace(0,0.3,21);

## Use the maximum coefficients, if available?
usemaxcoef=1;

#### Do NOT touch past here ####
if (!exist("octavedump.m","file"))
  error("This script requires an octavedump.m file to work.")
endif
source("octavedump.m");

for i = 1:columns(x)
  x(:,i) /= coef0(i);
endfor

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
    [w,iteration] = LassoActiveSet(x,y,t);
  endif

  worig = w;

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
  endif

  wrms = sqrt(sum((y - x * worig).^2));

  printf("%2.2d %.10f %.10f %.10f %.10f %.10f %d %d\n",it,t,...
         sum(abs(w)),sqrt(sum(w.^2)),max(abs(w)),...
         wrms,sum(sum(nterms)),iteration);
endfor

