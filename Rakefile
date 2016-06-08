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
@debug               = ENV.has_key?('DEBUG')
@user                = ENV['USER']
# basic compilers handled by the default configuration: ./config/default
@defaultCompilers    = %w[icpc icc clang clang++ gcc g++]
# default configure call
@defautConfigureCall = lambda {|cc| "./config/default CC=#{cc}"}
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
          puts line.strip.colorize(color: :green) if @debug and :out == key
          puts line.strip.colorize(color: :red)   if @debug and :err == key
        end
      end
    }

    # Don't exit until the external process is done
    external.join
  }
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
# execute remote command
def executeOnHost(command, builder)
  dbg(command)

  Net::SSH.start(builder.hostname,builder.username) {|ssh|
    out = ssh.exec!(command).strip
    dbg(out)
  }
end
#
# construct task from builder object
def builder2task(builder)
  baseTaskName   = "#{builder.host}#{builder.compiler.upcase}"
  syncTaskName   = "#{baseTaskName}_sync"
  configTaskName = "#{baseTaskName}_conf"
  buildTaskName  = "#{baseTaskName}_make"
  checkTaskName  = "#{baseTaskName}_check"

  desc "sync files for host: #{builder.host}, branch: #{getBranchName}"
  task syncTaskName.to_sym do |t|
    dbg("sync source  code for branch:" + getBranchName)
    doSync(builder)
  end

  desc "configure on host: %s, compiler %s, branch: %s" % [builder.host, builder.compiler, getBranchName]
  task configTaskName.to_sym do |t|
    dbg("call #{builder.configureCall}")
  end

  desc "build on host: %s, compiler %s, branch: %s" % [builder.host, builder.compiler, getBranchName]
  task buildTaskName.to_sym do |t|
    dbg("call 'make'")
  end

  desc "check on host: %s, compiler %s, branch: %s" % [builder.host, builder.compiler, getBranchName]
  task checkTaskName.to_sym do |t|
    dbg("call 'make check'")
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
# remote sync task {{{
# }}}
#
# remote configure task {{{
# }}}
#
# remote building task {{{
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
  dbg(getBranchName)
  dbg(syncFileList) if false
  dbg(executeOnHost("pwd",@userConfig["hosts"]["thunder4"]))
  dbg(executeOnHost("pwd",@userConfig["hosts"]["cygwin"]))
end
# }}}
