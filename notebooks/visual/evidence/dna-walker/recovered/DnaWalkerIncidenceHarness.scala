package hz.ricardo
package graph_rewriting
import moments._
object DnaWalkerIncidenceHarness {
  import DnaWalkerHarness.{N,E,NL,EL,Gr,Poly,R,G,graphBuilder,quote,rules,g2,g4}
  // Faithful incidence encoding: directed edge identity becomes a vertex.
  // This retains the permutations of two parallel, indistinguishable legs.
  def encode(g: Gr): Gr = {
    val old=g.nodes.toSeq.map(n => ("v:"+n,g(n).label))
    val en=g.edges.toSeq.map(e => ("edge:"+e.id,Some(if(g(e.source).label==Some("bimotor")) "edge_leg" else "edge_dna")))
    val h=G((old++en):_*)()
    for(e<-g.edges) {
      val a=("v:"+e.source) ~~> ("edge:"+e.id)
      val b=("edge:"+e.id) ~~> ("v:"+e.target)
      h.addEdges(Seq(a,b))
    }
    h
  }
  val ie=encode(g2); val ic=encode(g4)
  def encodeRule(r: R): R = {
    val l=encode(r.lhs); val h=encode(r.rhs)
    val ns=r.n.map {case(a,b)=>("v:"+a,"v:"+b)} ++ r.e.map {case(a,b)=>("edge:"+a.id,"edge:"+b.id)}
    val es=l.edges.toSeq.filter(e => ns.contains(e.source) && ns.contains(e.target))
      .map(e => e -> h.edges.find(f => f.source==ns(e.source) && f.target==ns(e.target) && h(f).label==l(e).label).get).toMap
    Rule(l,h,ns,es,r.rate)
  }
  val irules=rules.map(encodeRule)
  def incidenceClass(g: Gr): String = {
    val bs=g.nodes.filter(n=>g(n).label==Some("bimotor"))
    val cs=g.nodes.filter(n=>g(n).label==Some("chain"))
    if(bs.size!=1) return "zero"
    require(cs.size<5,"a pattern spanning a complete allowed DNA cycle needs another reducer")
    val es=g.nodes.filter(n=>g(n).label.exists(_.startsWith("edge_")))
    if(bs.size+cs.size+es.size!=g.nodes.size) throw new IllegalArgumentException("unknown node label")
    var legs=Vector.empty[(N,N)]
    var dna=Vector.empty[(N,N)]
    for(n<-es) {
      val in=g(n).incoming.toSeq; val out=g(n).outgoing.toSeq
      if(in.size>1 || out.size>1) return "zero"
      require(in.size==1 && out.size==1,"partial edge incidence needs another reducer: "+g)
      val a=in.head.source; val b=out.head.target
      if(g(n).label==Some("edge_leg")) {
        if(!bs(a) || !cs(b)) return "zero"
        legs=legs:+(a,b)
      } else {
        if(!cs(a) || !cs(b)) return "zero"
        dna=dna:+(a,b)
      }
    }
    if(g.edges.size!=2*es.size) return "zero"
    if(legs.size>2) return "zero"
    require(legs.size==2,"partial walker needs another reducer")
    if(dna.exists {case(a,b)=>a==b} || dna.distinct.size!=dna.size) return "zero"
    if(dna.groupBy(_._1).values.exists(_.size>1)||dna.groupBy(_._2).values.exists(_.size>1)) return "zero"
    if(dna.size>=cs.size) return "zero"
    var seen=Set(cs.head); var old=Set.empty[N]
    while(seen!=old) {old=seen; for((a,b)<-dna if seen(a)||seen(b)) seen=seen+a+b}
    require(seen.size==cs.size,"disconnected DNA needs another reducer")
    val a=legs(0)._2; val b=legs(1)._2
    if(a==b) "Craw" else if(dna.contains((a,b))||dna.contains((b,a))) "E" else "zero"
  }
  def incidenceInvariant(g: Gr): Option[Poly] = incidenceClass(g) match {
    case "zero" => Some(Pn.zero)
    case "E" => if(g iso ie) None else Some(Pn(Mn(ie)))
    case "Craw" => if(g iso ic) None else Some(Pn(Mn(ic)))
  }
  def solveIncidence(rs: List[R],seeds: List[Gr],budget: Int): String = {
    val eqs=generateMeanODEs(budget,rs,seeds,incidenceInvariant _)
    val odes=ODEPrinter(eqs).simplify.toVector
    def term(p: Poly): String=p.terms.flatMap {m=>
      require(m.denom.isEmpty&&m.numer.size==1,"linear pattern expression required")
      m.coef.terms.map {c=>require(c.denom.isEmpty,"polynomial rate required");
        "{\"coefficient\":"+c.coef+",\"rates\":"+c.numer.map(r=>quote(r.name)).mkString("[",",","]")+",\"coordinate\":"+quote(incidenceClass(m.numer.head))+"}"
      }
    }.mkString("[",",","]")
    val rows=odes.map(e=>"{\"coordinate\":"+quote(incidenceClass(e.lhs))+",\"terms\":"+term(e.rhs)+"}")
    "{\"ode_count\":"+odes.size+",\"equations\":"+rows.mkString("[",",","]")+"}"
  }
  def main(args: Array[String]): Unit = {
    val parts=for(r<-irules;q<-List(ie,ic)) yield
      "{\"rule\":"+quote(r.rate.name)+",\"coordinate\":"+quote(incidenceClass(q))+",\"output\":"+solveIncidence(List(r),List(q),1)+"}"
    println("DNA_INCIDENCE_JSON={\"single_query_run\":"+solveIncidence(irules,List(ie),12)+",\"individual_channels\":"+parts.mkString("[",",","]")+"}")
  }
}
