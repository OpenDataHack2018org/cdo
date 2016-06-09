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
if @userConfig.empty? then
  warn "No host information!!"
  exit(1)
end
# get setup from the environment
@debug               = true # == Rake.verbose ? true : false
@user                = ENV['USER']
# default configure call
@defautConfigureCall = lambda {|cc| "./config/default CC=#{cc}"}
# }}}

# helper methods {{{ ===========================================================
# general debugging output
def dbg(msg); pp msg if @debug; end

# return name of current branch
def getBranchName; `git branch`.split("\n").grep(/^\*/)[0].split[-1]; end

# wrapper for calling commands locally
#
# stdout is shown in debug mode only
# stderr is always shown
def executeLocal(cmd)
  system(cmd)
  return
  Open3.popen3(cmd) {|stdin, stdout, stderr, external|
    # read from stdout and stderr in parallel
    { :out => stdout, :err => stderr }.each {|key, stream|
      Thread.new do
        until (line = stream.gets).nil? do
          puts line.chomp                       if :out == key
          puts line.chomp.colorize(color: :red) if :err == key
        end
      end
    }
    # Don't exit until the external process is done
    external.join
  }
end

# execute remote command and collect
#
# stdout is shown in debug mode only
# stderr is always shown
def executeRemote(command, builder)
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
#
# execution wrapper
def execute(command, builder)
  # work in the target directory, ONLY
  command = ["test -f /etc/profile && source /etc/profile",
             "test -f .profile && source .profile",
             "test ! -d #{builder.targetDir} && mkdir -p #{builder.targetDir}",
             "cd #{builder.targetDir}",
             command].join(';')

  dbg(command)

  if builder.isLocal? then
    executeLocal(command)
  else
    executeRemote(command,builder)
  end
end
#
# synchronization for a given host
def doSync(builder)
  # make sure, that the remote target directory is present
  execute("mkdir -p #{builder.targetDir}", builder)

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
def builder2task(builder,useHostAsName=false,syncSource=true)
  baseTaskName    = useHostAsName ? builder.host : "#{builder.host}#{builder.compiler.upcase}"
  syncTaskName    = "#{baseTaskName}_sync"
  configTaskName  = "#{baseTaskName}_conf"
  buildTaskName   = "#{baseTaskName}_make"
  cleanTaskName   = "#{baseTaskName}_clean"
  checkTaskName   = "#{baseTaskName}_check"
  checkVTaskName  = "#{baseTaskName}_checkV"
  modlistTaskName = "#{baseTaskName}_mods"
  showLogTaskName = "#{baseTaskName}_showLog"

  if syncSource then
    #desc "sync files for host: #{builder.host}, branch: #{getBranchName}"
    task syncTaskName.to_sym do |t|
      dbg("sync source  code for branch:" + getBranchName)
      doSync(builder)
    end
  end

  #desc "configure on host: %s, compiler %s, branch: %s" % [builder.host, builder.compiler, getBranchName]
  task configTaskName.to_sym do |t|
    dbg("call #{builder.configureCall}")
    execute("#{builder.configureCall}",builder)
  end

  #desc "build on host: %s, compiler %s, branch: %s" % [builder.host, builder.compiler, getBranchName]
  task buildTaskName.to_sym do |t|
    execute("make -j4",builder)
  end

  #desc "check on host: %s, compiler %s, branch: %s" % [builder.host, builder.compiler, getBranchName]
  task checkTaskName.to_sym do |t|
    execute("make check",builder)
  end

  #desc "build on host: %s, compiler %s, branch: %s" % [builder.host, builder.compiler, getBranchName]
  task cleanTaskName.to_sym do |t|
    execute("make clean",builder)
  end

  #desc "check on host: %s, compiler %s, branch: %s" % [builder.host, builder.compiler, getBranchName]
  task checkVTaskName.to_sym do |t|
    execute("./src/cdo -V",builder)
  end

  # show remote config.log file
  task showLogTaskName.to_sym do |t|
    execute("cat config.log",builder)
  end

  # get the auto loaded modules on the target machine
  task modlistTaskName.to_sym do |t|
    execute("module list", builder)
  end

  desc "builder for host:#{builder.hostname}, CC=#{builder.compiler}"
  task baseTaskName.to_sym  => [syncSource ? syncTaskName : nil,
                                configTaskName,
                                buildTaskName,
                                checkTaskName].compact.map(&:to_sym)
end
# }}}
# constuct builders out of user configuration {{{ ==============================
Builder = Struct.new(:host,:hostname,:username,:compiler,:targetDir,:configureCall,:isLocal?)
# 1) construct builders from host configuration
#    this is what config/default should be able to handle
@userConfig["hosts"].each {|host,config|
  config['CC'].each {|cc|
    builder = Builder.new(host,
                          config["hostname"],
                          ('localhost' == config['hostname']) \
                              ? @user \
                              : ( config.has_key?('username') \
                                 ? config['username'] \
                                 : @userConfig["remoteUser"]),
                          cc,
                          [config["dir"],cc,getBranchName].join(File::SEPARATOR),
                          "./config/default CC=#{cc}",
                          config["hostname"] == 'localhost')
    builder2task(builder)
  }
}
# 2) construct builders from manual configuration
@userConfig["builders"].each {|builderName,config|
  builder = Builder.new(builderName,
                        @userConfig['hosts'][config["hostname"]]['hostname'],
                        ('localhost' == config['hostname'] \
                                  or 'localhost' == @userConfig['hosts'][config['hostname']]['hostname']) \
                            ? @user \
                            : ( config.has_key?('username') \
                               ? config['username'] \
                               : @userConfig["remoteUser"]),
                        '', # CC can be empty here, because it should be set by the given configureCall
                        [@userConfig['hosts'][config['hostname']]['dir'],builderName,getBranchName].join(File::SEPARATOR),
                        config['configureCall'],
                        ( 'localhost' == config['hostname'] \
                           or 'localhost' == @userConfig['hosts'][config['hostname']]['hostname'] ))
  builder2task(builder,true, config['sync'])

} if @userConfig.has_key?('builders')
# }}}
#

task :par do |t|
  # remove all tasks from the stack
  Rake.application.top_level_tasks.clear

  # create a task list from the command line
  ARGV.shift
  taskList = ARGV
  dbg(taskList)

  # execute tasks in parallel
  Parallel.map(taskList) {|t|
   sh "xterm -hold -e 'rake #{t}' "
  }
end
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
end
# }}}
