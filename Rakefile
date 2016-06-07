require 'json'
require 'net/ssh'
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
def syncFileList(file)
  file.write(`git ls-files`)
  dbg(File.open(file.path).readlines)
end
#
# synchronization for a given host
def doSync(builder)
  # remove eerything on the remote site first
  rsyncOpts = "--delete-excluded --delete"
  # keep old stuff on the remote site
  rsyncOpts = "-L"

  # collect the source files
  file = Tempfile.new("rsyncCdoTempfiles2Transfer")
  syncFileList(file)
  file.close

  # call rsync for a given host
  if builder.isLocal?
    syncCmd = "rsync #{rsyncOpts} -avz --files-from=#{file.path} . #{builder.targetDir}"
  else
    syncCmd = "rsync #{rsyncOpts} -avz --files-from=#{file.path}  -e '#{builder.command}' . #{builder.targetDir}"
  end
  dbg(syncCmd)
  call(syncCmd)
end
#
# execute remote command
def executeOnHost(command, config)
  dbg(command)

  hostname = config['hostname']
  username = 'localhost' == config["hostname"] \
    ? @user \
    : ( config.has_key?('username') \
       ? config['username'] \
       : @userConfig["remoteUser"] )

  Net::SSH.start(hostname,username) {|ssh|
    ssh.exec!(command).strip
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
Builder = Struct.new(:host,:hostname,:compiler,:targetDir,:configureCall,:configureOptions,:isLocal?)
@userConfig["hosts"].each {|host,config|
  @defaultCompilers.each {|cc|
    builder = Builder.new(host,
                          config["hostname"],
                          cc,
                          config["dir"],
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

        sh = ssh.shell.open

        pp sh.echo "#{config['dir']}"
        pp sh.mkdir(" -p #{config['dir']}")
        pp sh.cd("#{config['dir']}")
        puts sh.pwd.colorize(color: :red)
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
