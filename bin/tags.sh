#!/bin/bash
echo "Running tags.sh"
rm -f tags
for i in $(grep module *.F90 | sed '/end module/d' | cut -f1 -d.)
do
   mm=$(echo "$i" | cut -f1 -d\_)
   #tag=$(echo "$i" | cut -f2 -d\_)
   tag=$(echo "$i" | sed -e 's/m_//g' -e 's/mod_//g')
   echo "${tag}	${i}.F90	/^module"  >> tagsA
   echo "${i}	${i}.F90	/^module"  >> tagsA
done
export LC_COLLATE=C
cat tagsA  | sort -u --ignore-case > tags
rm -f tagsA
