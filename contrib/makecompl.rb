#!/usr/bin/env ruby
#== Synopsis
# Create auto completion files for different shells
#== Usage
#   makecompl.rb [-h] [-b]
#== Options
# -h    Show this help
# -b    Choose a different CDO binary for code generation 
#== Author
# Ralf Mueller, ralf.mueller@zmaw.de
#== RESTRICTIONS:
# TCSH: Completions is only performed for regular options and operators with prepended '-'
#== LICENSE: 
# CDO's License
#
require 'optparse'

options = {:bin => '../src/cdo'}
OptionParser.new do |opts|
  opts.banner = "Usage:  makecompl.rb [-h- [-b CMD]"
  opts.separator ""
  opts.separator "Options are ..."

opts.on("-h","--help","Show this help")  {puts opts}
opts.on("-b","--binary CMD","Choose a CDO binary different from #{options[:bin]}") {|cmd| options[:bin] = cmd}
end.parse!
#=============================================================================== 

outfile   = 'cdoCompletion'
cmd       = options[:bin] + ' 2>&1'
complCmds = { 
  :tcsh => ['set cdoCmpl = (\\','); complete cdo \'c/-/$cdoCmpl/\' \'n/*/f/\''],
  :zsh  => ['compctl -k "(',')" -f cdo'],
  :bash => ['complete -W "','" -f cdo']
}
help      = IO.popen(cmd).readlines.map {|l| l.chomp.lstrip}
options   = help.find_all {|item| /^-/.match(item)}.map {|o| o[0,2]}
operators = help[help.index("Operators:")+1].split
# (options + operators).each {|o| print o }; puts

[:bash, :zsh, :tcsh].each {|shell|
  completions = (:tcsh == shell ) ? options.map {|o| o[1..-1]}  + operators : options + operators.map {|o| "#{o} -#{o}"}
  File.open(outfile + '.' + shell.to_s,'w') {|f|
    f << complCmds[shell][0] << "\n"
    completions.each {|item| f << item << " \\" << "\n" }
    f << complCmds[shell][1] << "\n"
  }
}
