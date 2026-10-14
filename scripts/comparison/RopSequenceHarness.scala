package hz.ricardo
package graph_rewriting

import moments._

/** Positive directed-word graph input for the ungated ROP copolymer comparison.
  *
  * A state is one R -> word -> T path plus isolated F_A/F_B vertices.  This
  * harness intentionally supplies just the eight rules and the seven query
  * graphs below to FRAGGER.  It contains neither a suffix event partition nor
  * a target moment basis/equation.
  */
object RopSequenceHarness {
  type N = String
  type E = IdDiEdge[Int,N]
  type NL = String
  type EL = String
  type Gr = DiGraph[N,NL,E,EL]
  type Poly = Pn[N,NL,E,EL,DiGraph]
  val G = DiGraph.withType[N,NL,E,EL]
  implicit val graphBuilder = DiGraph.empty[N,NL,E,EL] _

  val Next = "next"
  val Root = "R"
  val Terminal = "T"

  def graph(nodes: Seq[(String,String)], arrows: Seq[(Int,String,String)]): Gr =
    G(nodes.map { case (id,label) => (id,Some(label)) }: _*)(
      // Fresh globally unique IDs are required by this upstream multigraph API:
      // local edge IDs from independently built patterns would otherwise clash
      // when FRAGGER forms unions.
      arrows.map { case (_,source,target) => (source ~~> target,Option.empty[String]) }: _*
    )

  def isolated(label: String, id: String = "v0"): Gr = graph(Seq(id -> label),Seq())

  /** A connected positive path, optionally with the unique R/T anchors. */
  def path(labels: Seq[String], prefixRoot: Boolean = false, suffixTerminal: Boolean = false): Gr = {
    val body = labels.zipWithIndex.map { case (label,index) => (s"w$index",label) }
    val nodes = (if (prefixRoot) Seq("r" -> Root) else Seq()) ++ body ++
      (if (suffixTerminal) Seq("t" -> Terminal) else Seq())
    val ids = (if (prefixRoot) Seq("r") else Seq()) ++ body.map(_._1) ++
      (if (suffixTerminal) Seq("t") else Seq())
    graph(nodes,ids.sliding(2).filter(_.size == 2).zipWithIndex.map {
      case (Seq(a,b),index) => (index,a,b)
    }.toSeq)
  }

  def seed(name: String, labels: String*): (String,Gr) = name -> path(labels)
  val seeds = Seq(
    seed("N_A","A"), seed("N_B","B"),
    seed("C_AA","A","A"), seed("C_AB","A","B"),
    seed("C_BA","B","A"), seed("C_BB","B","B"),
    seed("C_ABA","A","B","A")
  )

  def feed(symbol: String): Rule[N,NL,E,EL,DiGraph] =
    Rule(graph(Seq(),Seq()),isolated(s"F_$symbol"),Map(),Map(),Rate(s"beta_$symbol"))

  /** Preserve c and T, remove one isolated free pool vertex, and replace the
    * final c -> T edge by c -> symbol -> T.  R/A/B are concrete cases, so
    * the six rules are disjoint and exhaustive on the single-chain language.
    */
  def append(previous: String, symbol: String): Rule[N,NL,E,EL,DiGraph] = {
    val left = graph(Seq("c" -> previous,"t" -> Terminal,"f" -> s"F_$symbol"),Seq((0,"c","t")))
    val right = graph(Seq("c" -> previous,"x" -> symbol,"t" -> Terminal),Seq((1,"c","x"),(2,"x","t")))
    Rule(left,right,Map("c" -> "c","t" -> "t"),Map(),Rate(s"k_$symbol"))
  }

  val namedRules = Seq("feed_A" -> feed("A"),"feed_B" -> feed("B")) ++
    (for (symbol <- Seq("A","B"); previous <- Seq(Root,"A","B"))
      yield s"append_${previous}_$symbol" -> append(previous,symbol))
  val rules = namedRules.map(_._2)

  def quote(value: String): String = "\"" + value.flatMap {
    case '"' => "\\\""; case '\\' => "\\\\"; case '\b' => "\\b"; case '\f' => "\\f"
    case '\n' => "\\n"; case '\r' => "\\r"; case '\t' => "\\t"
    case character if character < ' ' => f"\\u${character.toInt}%04x"; case character => character.toString
  } + "\""
  def graphJson(value: Gr): String = {
    val nodes = value.nodes.toSeq.sorted.map(node =>
      "{\"id\":" + quote(node) + ",\"label\":" + quote(value(node).label.get) + "}")
    val edges = value.edges.toSeq.sortBy(_.id).map(edge =>
      "{\"id\":" + edge.id + ",\"source\":" + quote(edge.source) +
        ",\"target\":" + quote(edge.target) + ",\"label\":" + value(edge).label.map(quote).getOrElse("null") + "}"
    )
    "{\"nodes\":" + nodes.mkString("[",",","]") + ",\"edges\":" + edges.mkString("[",",","]") + "}"
  }
  def monomialJson(monomial: Mn[N,NL,E,EL,DiGraph]): String =
    monomial.graphs.map(graphJson).mkString("[",",","]")
  def polyJson(poly: Poly): String = poly.terms.flatMap { monomial =>
    require(monomial.denom.isEmpty,"no graph denominator in generated mean equation")
    monomial.coef.terms.map { rateTerm =>
      require(rateTerm.denom.isEmpty,"no rate denominator in generated mean equation")
      "{\"coefficient\":" + rateTerm.coef + ",\"rates\":" +
        rateTerm.numer.map(rate => quote(rate.name)).mkString("[",",","]") +
        ",\"graphs\":" + monomialJson(monomial) + "}"
    }
  }.mkString("[",",","]")

  def mapJson(values: Map[String,String]): String = values.toSeq.sortBy(_._1).map {
    case (left,right) => quote(left) + ":" + quote(right)
  }.mkString("{",",","}")
  def ruleJson(name: String,rule: Rule[N,NL,E,EL,DiGraph]): String =
    "{\"name\":" + quote(name) + ",\"rate\":" + quote(rule.rate.name) + ",\"lhs\":" + graphJson(rule.lhs) +
      ",\"rhs\":" + graphJson(rule.rhs) + ",\"preserved_nodes\":" + mapJson(rule.n) +
      ",\"preserved_edges\":" + rule.e.toSeq.sortBy(_._1.id).map { case (left,right) =>
        "{\"left\":" + left.id + ",\"right\":" + right.id + "}"
      }.mkString("[",",","]") + "}"

  def main(args: Array[String]): Unit = {
    if (args.headOption.contains("export")) {
      println("INPUT_JSON={\"edge_label_representation\":null,\"rules\":" + namedRules.map { case (name,rule) => ruleJson(name,rule) }.mkString("[",",","]") +
        ",\"seeds\":" + seeds.map { case (name,value) => "{\"name\":" + quote(name) + ",\"graph\":" + graphJson(value) + "}" }.mkString("[",",","]") + "}")
      return
    }
    val cap = args.headOption.map(_.toInt).getOrElse(12)
    val started = System.nanoTime()
    val equations = generateMeanODEs(cap,rules,seeds.map(_._2))
    val odes = ODEPrinter(equations).simplify.toVector
    val known = odes.map(_.lhs)
    val unresolved = odes.flatMap(_.rhs.graphs).filterNot(graph => known.exists(_ iso graph))
      .foldLeft(Vector.empty[Gr])((seen,graph) => if (seen.exists(_ iso graph)) seen else seen :+ graph)
    val rows = odes.map(eq => "{\"lhs\":" + graphJson(eq.lhs) + ",\"terms\":" + polyJson(eq.rhs) + "}")
    println("FRAGGER_JSON={\"ode_budget\":" + cap + ",\"seed_names\":" +
      seeds.map(value => quote(value._1)).mkString("[",",","]") + ",\"ode_count\":" + odes.size +
      ",\"algebraic_equations\":" + (equations.size-odes.size) + ",\"unresolved\":" +
      unresolved.map(graphJson).mkString("[",",","]") + ",\"seconds\":" +
      ((System.nanoTime()-started)/1e9) + ",\"equations\":" + rows.mkString("[",",","]") + "}")
  }
}
