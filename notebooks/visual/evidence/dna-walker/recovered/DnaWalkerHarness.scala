package hz.ricardo
package graph_rewriting

import moments._

/** JSON export of the published DNA-walker comparison.
  * The four rules and two transformers are those in upstream Bimotor.scala
  * at 0d9595c5d6e8918fc605f7d65778568a74f5e2d4. No target ODE is supplied.
  * Raw injective counts retain the two permutations of a compressed motor's
  * parallel leg edges. The Python comparison records physical normalization.
  */
object DnaWalkerHarness {
  type N = String
  type E = IdDiEdge[Int,N]
  type NL = String
  type EL = String
  type Gr = DiGraph[N,NL,E,EL]
  type Poly = Pn[N,NL,E,EL,DiGraph]
  type R = Rule[N,NL,E,EL,DiGraph]
  val G = DiGraph.withType[N,NL,E,EL]
  implicit val graphBuilder = DiGraph.empty[N,NL,E,EL] _
  val bc0 = "b" ~~> "c1"
  val bc1 = "b" ~~> "c1"
  val bc2 = "b" ~~> "c2"
  val bc3 = "b" ~~> "c2"
  val c1c2 = "c1" ~~> "c2"
  val g1 = G("b"->"bimotor","c1"->"chain","c2"->"chain")(bc0,bc1,c1c2)
  val g2 = G("b"->"bimotor","c1"->"chain","c2"->"chain")(bc1,bc2,c1c2)
  val g3 = G("b"->"bimotor","c1"->"chain","c2"->"chain")(bc2,bc3,c1c2)
  val g4 = G("b"->"bimotor","c1"->"chain")(bc0,bc1)
  val fe = Rule(g1,g2,Map("b"->"b","c1"->"c1","c2"->"c2"),Map(c1c2->c1c2,bc1->bc1),"kFE")
  val fc = Rule(g2,g3,Map("b"->"b","c1"->"c1","c2"->"c2"),Map(c1c2->c1c2,bc2->bc2),"kFC")
  val rules = List(fe,fc,fe.reversed("kBC"),fc.reversed("kBE"))

  def invariant(g: Gr): Option[Poly] =
    if (g.iso(g1) || g.iso(g3)) Some(Pn(Mn(g4))) else None
  def reachable(g: Gr): Option[Poly] = {
    val bs = g.nodes.filter(n => g(n).label == Some("bimotor"))
    val cs = g.nodes.filter(n => g(n).label == Some("chain"))
    if ((bs.size == 1) && (bs forall (g(_).outDegree == 2)) &&
        (cs.size <= 2) && (cs.toSeq.combinations(2) forall {
          case Seq(u,v) => ((g(u) edgesWith v).size == 1) })) None
    else Some(Pn.zero)
  }
  def quote(s: String): String = "\"" + s.replace("\\","\\\\").replace("\"","\\\"") + "\""
  def graphJson(g: Gr): String = {
    val ns=g.nodes.toSeq.sorted.map(n => "{\"id\":"+quote(n)+",\"label\":"+quote(g(n).label.get)+"}")
    val es=g.edges.toSeq.sortBy(_.id).map(e => "{\"id\":"+e.id+",\"source\":"+quote(e.source)+",\"target\":"+quote(e.target)+"}")
    "{\"nodes\":"+ns.mkString("[",",","]")+",\"edges\":"+es.mkString("[",",","]")+"}"
  }
  def ruleJson(r: R): String = {
    val ns=r.n.toSeq.sortBy(_._1).map {case(a,b) => "["+quote(a)+","+quote(b)+"]"}
    val es=r.e.toSeq.sortBy(_._1.id).map {case(a,b) => "["+a.id+","+b.id+"]"}
    "{\"name\":"+quote(r.rate.name)+",\"left\":"+graphJson(r.lhs)+",\"right\":"+graphJson(r.rhs)+
      ",\"preserved_nodes\":"+ns.mkString("[",",","]")+",\"preserved_edges\":"+es.mkString("[",",","]")+"}"
  }
  def polyJson(p: Poly): String = p.terms.flatMap { m =>
    require(m.denom.isEmpty && m.numer.size<=1,"linear graph expression expected")
    m.coef.terms.map { c =>
      require(c.denom.isEmpty,"polynomial rates expected")
      "{\"coefficient\":"+c.coef+",\"rates\":"+c.numer.map(r=>quote(r.name)).mkString("[",",","]")+
        ",\"graph\":"+(if(m.numer.isEmpty) "null" else graphJson(m.numer.head))+"}"
    }
  }.mkString("[",",","]")
  def solve(rs: List[R], seeds: List[Gr], budget: Int): String = {
    val eqs=generateMeanODEs(budget,rs,seeds,invariant _,reachable _)
    val odes=ODEPrinter(eqs).simplify.toVector
    val unresolved=odes.flatMap(_.rhs.graphs).filterNot(g=>odes.exists(e=>g iso e.lhs))
      .foldLeft(Vector.empty[Gr])((seen,g)=>if(seen.exists(_ iso g)) seen else seen :+ g)
    val rows=odes.map(e=>"{\"lhs\":"+graphJson(e.lhs)+",\"terms\":"+polyJson(e.rhs)+"}")
    "{\"ode_count\":"+odes.size+",\"equations\":"+rows.mkString("[",",","]")+
      ",\"unresolved\":"+unresolved.map(graphJson).mkString("[",",","]")+"}"
  }
  def main(args: Array[String]): Unit = {
    val pairs=for(r<-rules; q<-List(g2,g4)) yield
      "{\"rule\":"+quote(r.rate.name)+",\"query\":"+graphJson(q)+",\"output\":"+solve(List(r),List(q),1)+"}"
    println("DNA_WALKER_JSON={\"rules\":"+rules.map(ruleJson).mkString("[",",","]")+
      ",\"queries\":"+List(g1,g2,g3).map(graphJson).mkString("[",",","]")+
      ",\"extended_query\":"+graphJson(g2)+",\"compressed_query\":"+graphJson(g4)+
      ",\"original_query_run\":"+solve(rules,List(g1,g2,g3),12)+
      ",\"single_query_run\":"+solve(rules,List(g2),12)+
      ",\"individual_channels\":"+pairs.mkString("[",",","]")+"}")
  }
}
