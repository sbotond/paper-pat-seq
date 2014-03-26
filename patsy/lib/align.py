
import      utils       as      u
import      os
import      re
import      samtools    as      st

class Bowtie2:
    """ Bowtie2 wrapper class """
    def __init__(self, log, ref, rts, reads, index_opts={}, aln_opts={}, sink_err=True, has_index=False):
        self.prog       =       "bowtie2"
        self.log        =       log
        self.ref        =       ref
        self.rts        =       rts
        self.reads      =       reads
        self.index_opts =       index_opts
        self.aln_opts   =       aln_opts
        self.sink_err   =       sink_err
        self.has_index  =       has_index
        if not os.path.exists(ref):
            self.log.fatal(self.prog + ": Missing reference!")
        self.ref        =   os.path.abspath(ref)

    def index(self):
        """ Sort out indexing """
        # Symlink reference:
        new_ref         =   os.path.join(self.rts.base, os.path.basename(self.ref))
        os.symlink(self.ref,new_ref)
        # Register reference:
        self.rts.register(new_ref)

        old_ref         =   self.ref
        self.ref        =   new_ref
        ind_suf         =   [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]

        if self.has_index:
            # Symlink index:
            for sx in ind_suf:
                os.symlink(old_ref + sx, new_ref + sx)
        else:
            # Run indexing:
            cmd         = u.Cmd(self.log, prog="bowtie2-build",opts=self.index_opts,post_args=[self.ref, self.ref], sink_err=self.sink_err)
            cmd.comm()

        # Register index files:
        for suf in ind_suf:
            self.rts.register(self.ref + suf)
        pass

    @classmethod
    def mk_index(cls, ref, log, index_opts={}):
        """ Make a standalone index """
        cmd         = u.Cmd(log, prog="bowtie2-build",opts=index_opts,post_args=[ref, ref], sink_err=True)
        cmd.comm()

    def sam(self):
        """ Generate sam file """
        # Check input:
        if len(self.reads) == 0:
            self.log.fatal("No fastq files specified")

        (path, base)    = os.path.split(self.reads[0]) 
        pattern         = re.compile(".fq\d")
        base            = pattern.split(base)[0]

        # Generate SAM file:
        sam     = os.path.join(self.rts.base, base + ".sam")

        # Register sam file:
        self.rts.register(sam)

        # Construct post args:
        inp_flags   = ["-1", "-2"]
        inf         = [ ]
        for i in xrange(len(self.reads)):
            inf.append(inp_flags[i])            
            inf.append(self.reads[i])            

        # Add bwt base:
        tmp = [self.ref] 
        tmp.extend(inf)
        inf = tmp

        # Construct command object:
        cmd                 = u.Cmd(self.log, prog="bowtie2", opts=self.aln_opts, post_args=inf, outp=sam, cwd=self.rts.base, sink_err=self.sink_err)
        cmd.comm() 
        return sam


