require 'pp'

# Copyright (C) 2011-2012 Ralf Mueller, ralf.mueller@zmaw.de
# See COPYING file for copying and redistribution conditions.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# ==============================================================================
# CDO calling mechnism
module Cdo
  State = {
    :debug => false,
    :returnArray => false
  }
  @@CDO = ENV['CDO'].nil? ? '/usr/bin/cdo' : ENV['CDO']

  # Only operators with documentation are accessible vie the build-in help.
  # Other have to be added manually
  @@undocumentedOperators = %w[geopotheight pressure_fl pressure_hl]
  @@addOperators          = %w[boundaryLevels thicknessOfLevels]

  private
  def Cdo.call(cmd)
    if (State[:debug])
      puts '# DEBUG ====================================================================='
      puts cmd
      puts '# DEBUG ====================================================================='
      puts IO.popen(cmd).read
    else
      system(cmd + ' 1>/dev/null 2>&1 ')
    end
  end
  def Cdo.run(cmd,ofile=nil,options='')
    cmd = "#{@@CDO} -O #{options} #{cmd} "
    case ofile
    when $stdout
      cmd << " 2>/dev/null"
      return IO.popen(cmd).readlines.map {|l| l.chomp.strip}
    when nil
      ofile = Tempfile.new("Cdo.rb").path
    end
    cmd << "#{ofile}"
    call(cmd)
    if State[:returnArray]
      return NetCDF.open(ofile)
    else
      return ofile
    end
  end

  public
  def Cdo.Debug=(value)
    State[:debug] = value
  end
  def Cdo.Debug
    State[:debug]
  end
  def Cdo.returnArray=(value)
    if value
      begin
        require "numru/netcdf"
        include NumRu
        State[:returnArray] = true

      rescue LoadError
        warn "Could not load ruby's netcdf bindings. Please install it."
        raise
      end
    else
      State[:returnArray] = value
    end
  end
  def Cdo.returnArray
    State[:returnArray]
  end

  # test if @@CDO can be used
  def Cdo.checkCdo
    unless (File.exists?(@@CDO) and File.executable?(@@CDO))
      warn "Testing application #@@CDO is not available!"
      exit 1
    else
      puts "Using CDO: #@@CDO"
      puts IO.popen(@@CDO + " -V").readlines
    end
  end
  def Cdo.setCdo(cdo)
    puts "Will use #{cdo} instead of #@@CDO" if Cdo.Debug
    @@CDO = cdo
  end

  def Cdo.getOperators
    cmd       = @@CDO + ' 2>&1'
    help      = IO.popen(cmd).readlines.map {|l| l.chomp.lstrip}
    if 5 >= help.size
      warn "Operators could not get listed by running the CDO binary (#{@@CDO})"
      pp help if Cdo.Debug
      exit
    else
      help[help.index("Operators:")+1].split
    end
  end

  # Call an operator chain without checking opeartors
  def Cdo.chainCall(chain,*args)
    io = args.find {|a| a.class == Hash}
    args.delete_if {|a| a.class == Hash}

    chain   = chain.strip
    firstOp = chain
    firstOp = chain[0...[chain.index(','),chain.index(' ')].min] unless chain.index(',').nil?
    firstOp = firstOp[1..-1] if firstOp[0] == '-'
    if /(info|show|griddes)/.match(firstOp)
      Cdo.run(" #{chain} #{io[:in]} ",$stdout)
    else
      opts = args.empty? ? '' : ',' + args.reject {|a| a.class == Hash}.join(',')
      Cdo.run(" #{chain}#{opts} #{io[:in]} ",io[:out],io[:options])
    end
  end

  def Cdo.method_missing(sym, *args, &block)
    # args is expected to look like [opt1,...,optN,:in => iStream,:out => oStream] where
    # iStream could be another CDO call (timmax(selname(Temp,U,V,ifile.nc))
    puts "Operator #{sym.to_s} is called" if State[:debug]
    if getOperators.include?(sym.to_s) or @@undocumentedOperators.include?(sym.to_s)
      io = args.find {|a| a.class == Hash}
      args.delete_if {|a| a.class == Hash}
      if /(diff|info|show|griddes)/.match(sym)
        run(" -#{sym.to_s} #{io[:in]} ",$stdout)
      else
        opts = args.empty? ? '' : ',' + args.reject {|a| a.class == Hash}.join(',')
        run(" -#{sym.to_s}#{opts} #{io[:in]} ",io[:out],io[:options])
      end
    else
      warn "Operator #{sym.to_s} not found"
    end
  end

  #==================================================================
  # Addional operotors:
  #------------------------------------------------------------------
  def Cdo.boundaryLevels(args)
    ilevels         = Cdo.showlevel(:in => args[:in])[0].split.map(&:to_f)
    bound_levels    = Array.new(ilevels.size+1)
    bound_levels[0] = 0
    (1..ilevels.size).each {|i| 
      bound_levels[i] =bound_levels[i-1] + 2*(ilevels[i-1]-bound_levels[i-1])
    }
    bound_levels
  end

  def Cdo.thicknessOfLevels(args)
    bound_levels = Cdo.boundaryLevels(args)
    delta_levels    = []
    bound_levels.each_with_index {|v,i| 
      next if i == 0
      delta_levels << v - bound_levels[i-1]
    }
    delta_levels
  end
end

# Helper module for easy temp file handling
module MyTempfile
  require 'tempfile'
  @@_tempfiles           = []
  @@persistent_tempfiles = false
  def MyTempfile.setPersist(value)
    @@persistent_tempfiles = value
  end
  def MyTempfile.path
    unless @@persistent_tempfiles
      t = Tempfile.new(self.class.to_s)
      @@_tempfiles << t
      t.path
    else
      t = "_"+rand(10000000).to_s
      @@_tempfiles << t
      t
    end
  end
end
