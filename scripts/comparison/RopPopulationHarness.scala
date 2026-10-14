package hz.ricardo
package graph_rewriting

import moments._

/** Fresh FRAGGER worklist harness for the two small ROP population projections.
  *
  * This deliberately models only the declared common, edgeless population
  * language.  It is not an encoding of the molecular DPO graphs, nor does it
  * supply a target moment basis or any ODE.  ``baseline`` has Feed and
  * M -> N.  ``noedge`` has Feed, O <-> C, and O + M -> O + N.  All patterns
  * use injective embeddings and every displayed rule has per-match rate one.
  */
object RopPopulationHarness {
  type N = String
  type E = IdDiEdge[Int,N]
  type NL = String
  type EL = String
  type Gr = DiGraph[N,NL,E,EL]
  type Poly = Pn[N,NL,E,EL,DiGraph]
  val G = DiGraph.withType[N,NL,E,EL]
  implicit val graphBuilder = DiGraph.empty[N,NL,E,EL] _

  def population(labels: String*): Gr =
    G(labels.zipWithIndex.map { case (label,index) => (s"v$index",Some(label)) }: _*)()
  val empty = population()
  val incorporated = population("N")
  val monomer = population("M")
  val open = population("O")
  val closed = population("C")
  val openMonomer = population("O","M")
  val openIncorporated = population("O","N")

  def labels(graph: Gr): Seq[String] = graph.nodes.toSeq.map(node => graph(node).label.get).sorted

  /** The only supplied reduction for the NoEdge projection: O+C=1.
    * It is a state-language identity, not a list of desired output moments.
    */
  def singletonCatalyst(graph: Gr): Option[Poly] = {
    require(graph.edges.isEmpty,"the declared common language is edgeless")
    val found = labels(graph)
    val catalystCount = found.count(label => label == "O" || label == "C")
    if (catalystCount >= 2) Some(Pn.zero)
    else if (found.contains("C")) {
      val remaining = found.filterNot(_ == "C")
      Some(Pn(Mn(population(remaining: _*))) - Mn(population((remaining :+ "O"): _*)))
    } else None
  }

  def quote(value: String): String = "\"" + value.replace("\\","\\\\").replace("\"","\\\"") + "\""
  def graphJson(graph: Gr): String = labels(graph).map(quote).mkString("[",",","]")
  def polyJson(poly: Poly): String = poly.terms.flatMap { monomial =>
    require(monomial.denom.isEmpty,"no graph denominator in linear moment equations")
    require(monomial.numer.size <= 1,"do not silently factor joint moments")
    monomial.coef.terms.map { rational =>
      require(rational.denom.isEmpty,"no rate denominator")
      val graph = if (monomial.numer.isEmpty) "[]" else graphJson(monomial.numer.head)
      "{\"coefficient\":" + rational.coef + ",\"rates\":" +
        rational.numer.map(rate => quote(rate.name)).mkString("[",",","]") + ",\"graph\":" + graph + "}"
    }
  }.mkString("[",",","]")

  def run(variant: String): Unit = {
    val (rules, transformer) = variant match {
      case "baseline" => (
        List(
          Rule(empty,monomer,Map(),Map(),"beta_B"),
          Rule(monomer,incorporated,Map(),Map(),"k_B")
        ),
        None
      )
      case "noedge" => (
        List(
          Rule(empty,monomer,Map(),Map(),"beta_B"),
          Rule(openMonomer,openIncorporated,Map("v0"->"v0"),Map(),"k_B"),
          Rule(open,closed,Map(),Map(),"a"),
          Rule(closed,open,Map(),Map(),"b")
        ),
        Some(singletonCatalyst _)
      )
      case other => throw new IllegalArgumentException(s"unknown ROP population variant: $other")
    }
    // A common exploration allowance, deliberately larger than either
    // discovered closure.  It is not a requested basis dimension.
    val cap = 12
    val started = System.nanoTime()
    val equations = transformer match {
      case None => generateMeanODEs(cap,rules,List(incorporated))
      case Some(reduce) => generateMeanODEs(cap,rules,List(incorporated),reduce)
    }
    val odes = ODEPrinter(equations).simplify.toVector
    val unresolved = odes.flatMap(_.rhs.graphs).filterNot(graph => odes.exists(eq => graph iso eq.lhs))
      .foldLeft(Vector.empty[Gr])((seen,graph) => if (seen.exists(_ iso graph)) seen else seen :+ graph)
    val rows = odes.map(eq => "{\"lhs\":" + graphJson(eq.lhs) + ",\"terms\":" + polyJson(eq.rhs) + "}")
    println("FRAGGER_JSON={\"variant\":" + quote(variant) + ",\"seed\":[\"N\"],\"ode_budget\":" + cap +
      ",\"ode_count\":" + odes.size + ",\"algebraic_equations\":" + (equations.size-odes.size) +
      ",\"unresolved\":" + unresolved.map(graphJson).mkString("[",",","]") +
      ",\"seconds\":" + ((System.nanoTime()-started)/1e9) + ",\"equations\":" + rows.mkString("[",",","]") + "}")
  }

  def main(args: Array[String]): Unit = run(args.headOption.getOrElse("baseline"))
}
