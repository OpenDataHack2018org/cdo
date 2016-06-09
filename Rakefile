require 'json'
require 'net/ssh'
require 'open3'
require 'tempfile'
require 'pp'
require 'colorize'
require 'parallel'
require 'logger'

# configuration {{{ ============================================================
# load user setting if available
RC                   = "#{ENV['HOME']}/.rake.json"
@userConfig          = ( File.exist?(RC) ) ? JSON.load(File.open(RC)) : {}
# get setup from the environment
@debug               = true == Rake.verbose ? true : false
@user                = ENV['USER']
# basic compilers handled by the default configuration: ./config/default
@defaultCompilers    = %w[icpc icc clang clang++ gcc g++]
# default configure call
@defautConfigureCall = lambda {|cc| "./config/default CC=#{cc}"}
@srcDir              = File.expand_path(File.dirname(__FILE__))
# }}}

# helper methods {{{ ===========================================================
# general debugging output
def dbg(msg); pp msg if @debug; end
#
# collect file list for sync
def getBranchName
  `git branch`.split("\n").grep(/^\*/)[0].split[-1]
end
#
# get files for sync
def executeLocal(cmd)
  Open3.popen3(cmd) {|stdin, stdout, stderr, external|
    # read from stdout and stderr in parallel
    { :out => stdout, :err => stderr }.each {|key, stream|
      Thread.new do
        until (line = stream.gets).nil? do
          puts line.chomp                       if :out == key and @debug
          puts line.chomp.colorize(color: :red) if :err == key
        end
      end
    }

    # Don't exit until the external process is done
    external.join
  }
end
#
# execute remote command
def executeOnHost(command, builder)
  dbg(command)

  command = ["[[ -f /etc/process ]] && source /etc/profile",
             "[[ -f .profile ]] && source .profile",
             "cd #{builder.targetDir}",
             command].join(';')
  if builder.isLocal? then
    executeLocal(command)
  else
    Net::SSH.start(builder.hostname,builder.username) do |ssh|
      stdout_data = ""
      stderr_data = ""
      exit_code = nil
      exit_signal = nil
      ssh.open_channel do |channel|
        channel.exec(command) do |ch, success|
          unless success
            raise "FAILED: couldn't execute command #{command}"
          end
          channel.on_data do |ch, data|
            stdout_data += data
            $stdout.write(data) if @debug
          end

          channel.on_extended_data do |ch, type, data|
            stderr_data += data
            $stderr.write(data)# if @debug
          end

          channel.on_request("exit-status") do |ch, data|
            exit_code = data.read_long
          end

          channel.on_request("exit-signal") do |ch, data|
            exit_signal = data.read_long
          end
        end
      end
      ssh.loop
    end
  end
end
#
# synchronization for a given host
def doSync(builder)
  # make sure, that the remote target directory is present
  executeOnHost("mkdir -p #{builder.targetDir}", builder)

  # basic rsync options
  # * remove everything on the remote site first
  rsyncOpts = "--delete-excluded --delete"
  # * keep old stuff on the remote site
  rsyncOpts = "-L"

  # collect the source files
  file = Tempfile.new("rsyncCdoTempfiles2Transfer")
  begin
    file.write(`git ls-files`)
    file.write(`git submodule foreach 'git ls-files | sed "s|^|$path/|"'`)

    # call rsync for a given host
    if builder.isLocal?
      syncCmd = "rsync #{rsyncOpts} -avz --files-from=#{file.path} . #{builder.targetDir}"
    else
      syncCmd = "rsync #{rsyncOpts} -avz --files-from=#{file.path}  -e ssh . #{builder.username}@#{builder.hostname}:#{builder.targetDir}"
    end
    dbg(syncCmd)
    executeLocal(syncCmd)
  ensure
    file.close
    file.unlink
  end
end
#
# construct task from builder object
def builder2task(builder)
  baseTaskName   = "#{builder.host}#{builder.compiler.upcase}"
  syncTaskName   = "#{baseTaskName}_sync"
  configTaskName = "#{baseTaskName}_conf"
  buildTaskName  = "#{baseTaskName}_make"
  cleanTaskName  = "#{baseTaskName}_clean"
  checkTaskName  = "#{baseTaskName}_check"
  checkVTaskName = "#{baseTaskName}_checkV"

  desc "sync files for host: #{builder.host}, branch: #{getBranchName}"
  task syncTaskName.to_sym do |t|
    dbg("sync source  code for branch:" + getBranchName)
    doSync(builder)
  end

  desc "configure on host: %s, compiler %s, branch: %s" % [builder.host, builder.compiler, getBranchName]
  task configTaskName.to_sym do |t|
    dbg("call #{builder.configureCall}")
    executeOnHost("cd #{builder.targetDir}; #{builder.configureCall}",builder)
  end

  desc "build on host: %s, compiler %s, branch: %s" % [builder.host, builder.compiler, getBranchName]
  task buildTaskName.to_sym do |t|
    executeOnHost("cd #{builder.targetDir}; make -j4",builder)
  end

  desc "build on host: %s, compiler %s, branch: %s" % [builder.host, builder.compiler, getBranchName]
  task cleanTaskName.to_sym do |t|
    executeOnHost("cd #{builder.targetDir}; make clean",builder)
  end

  desc "check on host: %s, compiler %s, branch: %s" % [builder.host, builder.compiler, getBranchName]
  task checkTaskName.to_sym do |t|
    executeOnHost("cd #{builder.targetDir}; make check",builder)
  end

  desc "check on host: %s, compiler %s, branch: %s" % [builder.host, builder.compiler, getBranchName]
  task checkVTaskName.to_sym do |t|
    executeOnHost("cd #{builder.targetDir}; ./src/cdo -V",builder)
  end

  desc "builder for host:#{builder.hostname}, CC=#{builder.compiler}"
  task baseTaskName.to_sym  => [syncTaskName, configTaskName, buildTaskName, checkTaskName].map(&:to_sym) do |t|
    pp builder.to_h
  end
end
# }}}
# constuct builders out of user configuration {{{ ==============================
Builder = Struct.new(:host,:hostname,:username,:compiler,:targetDir,:configureCall,:configureOptions,:isLocal?)
@userConfig["hosts"].each {|host,config|
  compilers = config.has_key?('CC') ? config['CC'] : @defaultCompilers
  compilers.each {|cc|
    builder = Builder.new(host,
                          config["hostname"],
                          ('localhost' == config['hostname']) \
                              ? @user \
                              : ( config.has_key?('username') \
                                 ? config['username'] \
                                 : @userConfig["remoteUser"]),
                          cc,
                          [config["dir"],cc,getBranchName].join(File::SEPARATOR),
                          @defautConfigureCall[cc],
                          "",
                          config["hostname"] == 'localhost')
    builder2task(builder)
  }
}
# }}}
#
# check connections {{{
desc "check available connections"
task :checkConnections do |t|
  pp Parallel.map(@userConfig["hosts"]) {|host, config|
    hostname = config['hostname']
    username = 'localhost' == config["hostname"] \
      ? @user \
      : ( config.has_key?('username') \
         ? config['username'] \
         : @userConfig["remoteUser"] )

    if 'localhost' == config["hostname"] then
      config["hostname"]
    else
      Net::SSH.start(hostname,username) do |ssh|
        hn = ssh.exec!("hostname -f").strip
        us = ssh.exec!("echo $USER").strip
        [hn,us]
      end
    end
  }
end
# }}}
#
# check internals {{{
desc "check some internals"
task :checkInterals do
  dbg(@srcDir)
  dbg(getBranchName)
  dbg(syncFileList) if false
# dbg(executeOnHost("pwd",@userConfig["hosts"]["thunder4"]))
# dbg(executeOnHost("pwd",@userConfig["hosts"]["cygwin"]))
end
# }}}
