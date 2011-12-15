$:.unshift File.join(File.dirname(__FILE__),"..","lib")
require 'test/unit'
require 'cdo'

class TestJobQueue < Test::Unit::TestCase

  DEFAULT_CDO_PATH = '/usr/bin/cdo'
  def setup
    if ENV['CDO'].nil?
      if File.exists?(DEFAULT_CDO_PATH)
        Cdo.setCdo(DEFAULT_CDO_PATH)
      else
        stop(DEFAULT_CDO_PATH)
      end
    else
      # Check user given path
      unless File.exists?(ENV['CDO'])
        stop(ENV['CDO'])
      else
        Cdo.setCdo(ENV['CDO'])
      end
    end
  end
  def stop(path)
    warn "Could not find CDO binary (#{path})! Abort tests"
    exit
  end

  def test_getOperators
    %w[for random stdatm info showlevel sinfo remap mask topo thicknessOfLevels].each {|op|
      assert(Cdo.getOperators.include?(op),"Operator '#{op}' not found") unless op == "thicknessOfLevels"
      assert(Cdo.respond_to?(op),"Operator '#{op}' not found") if op == "thicknessOfLevels"
    }
  end
  def test_info
    levels = Cdo.showlevel(:in => "-stdatm,0")
    assert_equal([0,0].map(&:to_s),levels)

    info = Cdo.sinfo(:in => "-stdatm,0")
    assert_equal("File format: GRIB",info[0])

  end
  def test_args
    #Cdo.Debug = true
    #MyTempfile.setPersist(true)
    ofile0 = MyTempfile.path
    ofile1 = MyTempfile.path
    ofile2 = MyTempfile.path
    ofile3 = MyTempfile.path
    Cdo.stdatm(0,20,40,80,200,230,400,600,1100,:out => ofile0)
    Cdo.intlevel(0,10,50,100,500,1000,  :in => ofile0,:out => ofile1)
    Cdo.intlevel([0,10,50,100,500,1000],:in => ofile0,:out => ofile2)
    Cdo.sub(:in => [ofile1,ofile2].join(' '),:out => ofile3)
    info = Cdo.infon(:in => ofile3)
    (1...info.size).each {|i| assert_equal(0.0,info[i].split[-1].to_f)}
  end
  def test_operator_options
    ofile = MyTempfile.path
    targetLevels = [0,10,50,100,200,400,1000]
    Cdo.stdatm(targetLevels,:out => ofile)
    levels = Cdo.showlevel(:in => ofile)
    [0,1].each {|i| assert_equal(targetLevels.map(&:to_s),levels[i].split)}
  end
  def test_CDO_options
    names = Cdo.showname(:in => "-stdatm,0",:options => "-f nc")
    assert_equal(["P T"],names)

    ofile = MyTempfile.path
    Cdo.topo(:out => ofile,:options => "-z szip")
    assert_equal(["GRIB SZIP"],Cdo.showformat(:in => ofile))
  end
  def test_chain
    ofile     = MyTempfile.path
    #Cdo.Debug = true
    Cdo.chainCall("-setname,veloc -copy",:in => "-random,r1x1",:out => ofile,:options => "-f nc")
    assert_equal(["veloc"],Cdo.showname(:in => ofile))
  end

  def test_diff
    diffv = Cdo.diffv(:in => "-random,r1x1 -random,r1x1")
    assert_equal(diffv[1].split(' ')[4],"random")
    assert_equal(diffv[1].split(' ')[-1],"0.53060")
    diff  = Cdo.diff(:in => "-random,r1x1 -random,r1x1")
    assert_equal(diff[1].split(' ')[-1],"0.53060")
  end
end

#  # Calling simple operators
#  #
#  # merge:
#  #   let files be an erray of valid filenames and ofile is a string
#  Cdo.merge(:in => outvars.join(" "),:out => ofile)
#  #   or with multiple arrays:
#  Cdo.merge(:in => [ifiles0,ifiles1].flatten.join(' '),:out => ofile)
#  # selname:
#  #   lets grep out some variables from ifile:
#  ["T","U","V"].each {|varname|
#    varfile = varname+".nc"
#    Cdo.selname(varname,:in => ifile,:out => varfile)
#  }
#  #   a threaded version of this could look like:
#  ths = []
#  ["T","U","V"].each {|outvar|
#    ths << Thread.new(outvar) {|ovar|
#      varfile = varname+".nc"
#      Cdo.selname(varname,:in => ifile,:out => varfile)
#    }
#  }
#  ths.each {|th| th.join}
#  # another example with sub:
#  Cdo.sub(:in => [oldfile,newfile].join(' '), :out => diff)
#  
#  # It is possible too use the 'send' method
#  operator  = /grb/.match(File.extname(ifile)) ? :showcode : :showname
#  inputVars = Cdo.send(operator,:in => ifile)
#  # show and info operators are writing to stdout. cdo.rb tries to collects this into arrays
#  #
#  # Same stuff with other operators:
#  operator = case var
#             when Fixnum then 'selcode'
#             when String then 'selname'
#             else
#               warn "Wrong usage of variable identifier for '#{var}' (class #{var.class})!"
#             end
#  Cdo.send(operator,var,:in => @ifile, :out => varfile)
#  
#  # For chaining operators, there is a special method:
#  Cdo.chainCall("-setname,veloc -copy",:in => ifile,:out => ofile,:options => "-f nc")
#  #   another example with 3 operators and a different hash syntax
#  C_R          = 287.05
#  Cdo.chainCall("setname,#{rho} -divc,#{C_R} -div",in: [pressureFile,temperatureFile].join(' '), out: densityFile)
#  
#  # Pass an array for operators with multiple options:
#  #   Perform conservative remapping with pregenerated weights
#  Cdo.remap([gridfile,weightfile],:in => copyfile,:out => outfile)
#  #   Create vertical height levels out of hybrid model levels
#  Cdo.ml2hl([0,20,50,100,200,400,800,1200].join(','),:in => hybridlayerfile, :out => reallayerfile)
#  # or use multiple arguments directly
#  Cdo.remapeta(vctfile,orofile,:in => ifile,:out => hybridlayerfile)
#  
#  # the powerfull expr operator:
#  # taken from the tutorial in https://code.zmaw.de/projects/cdo/wiki/Tutorial#The-_expr_-Operator
#  SCALEHEIGHT  = 10000.0
#  C_EARTH_GRAV = 9.80665
#  # function for later computation of hydrostatic atmosphere pressure
#  PRES_EXPR    = lambda {|height| "101325.0*exp((-1)*(1.602769777072154)*log((exp(#{height}/#{SCALEHEIGHT})*213.15+75.0)/288.15))"}
#  TEMP_EXPR    = lambda {|height| "213.0+75.0*exp(-#{height}/#{SCALEHEIGHT})"}
#  
#  # Create Pressure and Temperature out of a height field 'geopotheight' from ifile
#  Cdo.expr("'p=#{PRES_EXPR['geopotheight']}'", :in => ifile, :out => presFile)
#  Cdo.expr("'t=#{TEMP_EXPR['geopotheight']}'", :in => ifile, :out => tempFile)
#  
#  
#  # TIPS: I often work with temporary files and for getting rid of handling them manually the MyTempfile module can be used:
#  #       Simply include the following methods into you scripts and use tfile for any temporary variable
#  def tfile
#    MyTempfile.path
#  end
#  # As an example, the computation of simple atmospherric density could look like
#  presFile, tempFile = tfile, tfile
#  Cdo.expr("'p=#{PRES_EXPR['geopotheight']}'", :in => ifile, :out => presFile)
#  Cdo.expr("'t=#{TEMP_EXPR['geopotheight']}'", :in => ifile, :out => tempFile)
#  Cdo.chainCall("setname,#{rho} -divc,#{C_R} -div",in: [presFile,tempFile].join(' '), out: densityFile)
#  
#  # For debugging, it is helpfull, to avoid the automatic cleanup at the end of the scripts:
#  MyTempfile.setPersist(true)
#  # creates randomly names files. Switch on debugging with 
#  Cdo.Debug = true
