package hz.ricardo
package graph_rewriting

import moments._

object ABScalingHarness {
  type N = String
  type E = IdDiEdge[Int,N]
  type NL = String
  type EL = String
  type G0 = DiGraph[N,NL,E,EL]
  type GG[X,Y,Z <: DiEdgeLike[X],W] = DiGraph[X,Y,Z,W]

  val G = DiGraph.withType[N,NL,E,EL]
  implicit val graphBuilder = DiGraph.empty[N,NL,E,EL] _

  def graph(nodes: Seq[(String,String)], arrows: Seq[(String,String)]): G0 =
    G(nodes.map { case (id,label) => (id,Some(label)) }: _*)(
      arrows.map { case (source,target) =>
        (source ~~> target, Option.empty[String])
      }: _*
    )

  def wordGraph(word: String): G0 = {
    val nodes = word.zipWithIndex.map { case (symbol,index) =>
      ("w" + index, symbol.toString)
    }
    val arrows = (0 until math.max(0, word.length - 1)).map { index =>
      ("w" + index, "w" + (index + 1))
    }
    graph(nodes, arrows)
  }

  def append(previous: String, symbol: String): Rule[N,NL,E,EL,GG] = {
    val left = graph(
      Seq("p" -> previous, "t" -> "T"),
      Seq("p" -> "t")
    )
    val right = graph(
      Seq("p" -> previous, "x" -> symbol, "t" -> "T"),
      Seq("p" -> "x", "x" -> "t")
    )
    Rule(left, right, Rate("append_" + previous + "_" + symbol))
  }

  val namedRules: Vector[(String,Rule[N,NL,E,EL,GG])] =
    (for {
      symbol <- Vector("A","B")
      previous <- Vector("R","A","B")
    } yield ("append_" + previous + "_" + symbol) -> append(previous,symbol))

  case class RawResult(total: Int, perOrder: Vector[(String,String,Int)])

  def enumerate(observable: G0): RawResult = {
    val rows = Vector.newBuilder[(String,String,Int)]
    var total = 0
    namedRules.foreach { case (name,rule) =>
      val left = observable.leftUnions(rule).size
      rows += ((name,"observable-left",left))
      total += left
      val right = observable.rightUnions(rule).size
      rows += ((name,"observable-right",right))
      total += right
    }
    RawResult(total,rows.result())
  }

  def median(values: Vector[Double]): Double = {
    val sorted = values.sorted
    if (sorted.length % 2 == 1) sorted(sorted.length / 2)
    else (sorted(sorted.length / 2 - 1) + sorted(sorted.length / 2)) / 2.0
  }

  val Q = 34.toChar.toString

  def jsonString(value: String): String = Q + value + Q

  def orderJson(rows: Vector[(String,String,Int)]): String =
    rows.map { case (rule,order,count) =>
      "{" +
      Q + "rule" + Q + ":" + jsonString(rule) + "," +
      Q + "order" + Q + ":" + jsonString(order) + "," +
      Q + "count" + Q + ":" + count +
      "}"
    }.mkString("[",",","]")

  def benchmark(family: String, n: Int, warmups: Int, repetitions: Int): Unit = {
    val word =
      if (family == "same") "A" * n
      else ("AB" * ((n + 1) / 2)).take(n)
    val observable = wordGraph(word)

    var reference: RawResult = null
    var warm = 0
    while (warm < warmups) {
      reference = enumerate(observable)
      warm += 1
    }

    val timings = Vector.newBuilder[Double]
    var result: RawResult = null
    var rep = 0
    while (rep < repetitions) {
      val start = System.nanoTime()
      result = enumerate(observable)
      val stop = System.nanoTime()
      require(result == reference, "non-deterministic raw composition result")
      timings += ((stop - start).toDouble / 1e9)
      rep += 1
    }
    val times = timings.result()

    val json =
      "{" +
      Q + "tool" + Q + ":" + jsonString("fragger-core") + "," +
      Q + "family" + Q + ":" + jsonString(family) + "," +
      Q + "word" + Q + ":" + jsonString(word) + "," +
      Q + "n" + Q + ":" + n + "," +
      Q + "raw_composites" + Q + ":" + result.total + "," +
      Q + "per_order" + Q + ":" + orderJson(result.perOrder) + "," +
      Q + "times_seconds" + Q + ":[" + times.mkString(",") + "]," +
      Q + "median_seconds" + Q + ":" + median(times) + "," +
      Q + "min_seconds" + Q + ":" + times.min + "," +
      Q + "max_seconds" + Q + ":" + times.max +
      "}"

    println("RESULT=" + json)
  }

  def main(args: Array[String]): Unit = {
    if (args.length < 2) {
      System.err.println("usage: ABScalingHarness FAMILY N [WARMUPS] [REPETITIONS]")
      sys.exit(2)
    }
    val family = args(0)
    require(family == "alt" || family == "same", "family must be alt or same")
    val n = args(1).toInt
    val warmups = if (args.length > 2) args(2).toInt else 2
    val repetitions = if (args.length > 3) args(3).toInt else 5

    val global = wordGraph("ABABA")
    var i = 0
    while (i < 4) {
      enumerate(global)
      i += 1
    }

    benchmark(family,n,warmups,repetitions)
  }
}
