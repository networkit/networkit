#!/bin/bash
RETURNCODE=0

echo -e "We will execute a number of check intended to improve code quality."
echo -e "Typically violations can be automatically fixed by executing:"
echo -e "./check_code.sh -w"
echo -e "\e[31m\e[1mWhile we're trying to be careful ALWAYS keep a copy of the original around!\e[39m\e[0m"
echo

for check in CppIndentation.py CppClangFormat.py CppIncludeGuards.py; do
  echo -e "\e[32m=== Execute $check ===\e[39m"
  if extrafiles/tooling/$check $@ ; then
    echo -e "\e[32mOK.\e[39m"
  else
    echo -e "\e[31mERROR.\e[39m"
    RETURNCODE=1
  fi
  echo
done

exit $RETURNCODE