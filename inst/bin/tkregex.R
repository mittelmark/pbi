library(tcltk)
.ttop=tktoplevel()
tkwm.title(.ttop,"Detlef's Regexer")
.thits=tclVar("")
.tregex=tclVar("")
.thpros=tclVar("")
.pbi.tmpros <- function () {
    pattern=tclvalue(tkget(tkps))
    pattern=gsub(" --.+","",pattern)
    pattern=gsub("^<","^",pattern)
    pattern=gsub(">$","$",pattern)
    pattern=gsub("\\{([^\\}]+)\\}","[^\\1]",pattern)
    pattern=gsub("\\((.+?)\\)","{\\1}",pattern)
    pattern=gsub("-","",pattern)
    pattern=gsub("x",".",pattern)
    tclvalue(.tregex)=pattern
    .pbi.tmfunc()
    #tkdelete(tke,'1','end')
      #tkinsert(tke,'1',pattern)
  }
  .pbi.tmfunc <<- function () {
      text=tclvalue(tkget(tke))
      if (text == "") {
          return()
      }
      tktag.delete(tkt,"found",'1.0','end')
      tktag.configure(tkt,"found",foreground="#aa3333")
                      c=tclVar("")
                      lidx=strsplit(tclvalue(tksearch(tkt,'-all','-regex',text,"1.0"))," ")[[1]]
                      tclvalue(.thits)=length(lidx)
                      for (idx in lidx) {
                          x=tclvalue(tksearch(tkt,'-regex',count=c,text,idx))
                          col=as.numeric(gsub("(.+)\\.([0-9]+)","\\2",idx))
                          lin=as.numeric(gsub("(.+)\\.([0-9]+)","\\1",idx))
                          end=paste(lin,col+as.numeric(tclvalue(c)),sep=".")
                          tktag.add(tkt,"found",idx,end)
                          
                      }
                  }
                  ftop=tkframe(.ttop)
                  tkl=tklabel(ftop,text="Prosite Expression:\npress Enter ...)")
                  tkps=tkentry(ftop,width=30)
    tkinsert(tkps,'1',"<M[AC]-x(0,1)-V-x(4)-{ED}-x(4,9)-M>")
    tkbind(tkps,'<Return>',.pbi.tmpros)
    tkgrid(tkl,padx=5,pady=5,row=0,column=0)
    tkgrid(tkps,padx=5,pady=5,row=0,column=1,sticky="nsew")
    
    tkl=tklabel(ftop,text="Regular Expression:\n(press Enter ...)")
    tke=tkentry(ftop,width=30,textvariable=.tregex)
    tkbind(tke,'<Return>',.pbi.tmfunc)
    tkgrid(tkl,padx=5,pady=5,row=1,column=0)
    tkgrid(tke,padx=5,pady=5,row=1,column=1,sticky='nsew')
    tkgrid(tklabel(ftop,text="Hits: "),row=0,column=2,padx=5,pady=5,rowspan=2,sticky="nsew")
    tkgrid(tklabel(ftop,text=".....",width=10,textvariable=.thits),row=0,
           column=3,padx=5,pady=5,rowspan=2)
    tkgrid(tkbutton(ftop,text="Close Application",
                    command=function() { tkdestroy(.ttop); }),
           row=0,column=4,rowspan=2,sticky='nsew',padx=20,pady=10)
    tkpack(ftop,side='top',fill='x',expand=FALSE)
    tkt=tktext(.ttop,border=5,relief='flat')
    tktag.configure(tkt,"found",foreground="#aa3333")
    tkpack(tkt,side='top',fill='both',expand=TRUE)
    tkinsert(tkt,"end",">ID1\nMAVFGTRTGHYRTGVM\n")
    tkinsert(tkt,"end",">ID2\nAFGTRTGHZKLMO\n")
tkwait.window(.ttop)
