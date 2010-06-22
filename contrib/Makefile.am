EXTRA_DIST = makecompl.rb cdoCompletion.bash cdoCompletion.tcsh cdoCompletion.zsh

all-am:
	if test -x /usr/bin/ruby -o -x /usr/local/bin/ruby ; then \
	  ./makecompl.rb; \
	fi
