#! /bin/bash

extract_block.awk empty.log 1 1 | grep '^ *|' | awk '{print $6}' > ref.dat
extract_block.awk empty.log 1 1 | grep '^ *|' | awk '{print $8}' > empty.dat

for k in {1..20} ; do
    extract_block.awk h.log $k $k | grep '^ *|' | awk '{print $8}'  > h_l_${k}.dat
    extract_block.awk h.log $((k+20)) $((k+20)) | grep '^ *|' | awk '{print $8}'  > h_s_${k}.dat
done

for i in {c,n,o,f,p,s,cl} ; do
    for k in {1..20} ; do
	extract_block.awk ${i}.log $k $k | grep '^ *|' | awk '{print $8}'  > ${i}_l_${k}.dat
	extract_block.awk ${i}.log $((k+20)) $((k+20)) | grep '^ *|' | awk '{print $8}'  > ${i}_s_${k}.dat
	extract_block.awk ${i}.log $((k+40)) $((k+40)) | grep '^ *|' | awk '{print $8}'  > ${i}_p_${k}.dat
    done
done

