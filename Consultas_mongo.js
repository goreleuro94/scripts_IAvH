/*Consultar los modelos de un taxon especifico y un status especifico
Estan agrupados los campos menos el _id para que no se dupliquen ya que los taxones que tienen sinonimias tienen el mismo taxonID por lo que al momento de hacer la consulta trae los dos
nombres como registros separados
*/

db.models.aggregate([
    {
        $lookup: {
            from: "species",
            localField: "taxID",
            foreignField: "taxID",
            as: "especie"
        }
    },
    {
        $unwind: "$especie"
    },
    {
        $match: {
            "modelStatus": "Valid",
            "especie.class": "Amphibia"
        }
    },
    {
        $group:
        { _id:"$_id",
          taxID: { $first: "$taxID" },
          modelID: { $first: "$modelID" },
          modelStatus: { $first: "$modelStatus" },
          kingdom: { $first: "$especie.kingdom" },
          order: { $first: "$especie.order" },
          class: { $first: "$especie.class" },
          family: { $first: "$especie.family" },
          genus: { $first: "$especie.genus" },      
          acceptedNameUsage: { $first: "$acceptedNameUsage" }  
        }
        }
])

/*Buscar por modelID en models*/
db.getCollection("models").find({modelID:"HER-426"})

/*Buscar por taxId en species*/
db.getCollection("species").find({taxID:7265})

/*Buscar registros por un taxID*/
db.getCollection("records").find({taxID:7265})

/*Buscar registros por varios taxID*/
db.getCollection("records").find({taxID:{ $in: [218,8139,241,8151,8156,5239,467,470,552,5255,5275,7093,7094,865,979,5375,1041,8189,1110,8195,8478,1187,1188,1298,7098,1410,5565,1600
    ,8479,8233,1811,2011,8248,2016,8247,2101,2102,7132,8258,8259,2318,8282,6993,2627,2692,8293,2782,3051,3062,6181,567,7115,3238,7118,8474,3288,8337,6272,8338,6328,3668,3690,8358,
    7123,4063,4065,4076,8382,153,154,8389,8390,4320,8395,8480,7126,4626,4629,6995,8453,8457] } }).itcount()

/*Buscar registros por acceptedNameUsage*/
db.getCollection('species').aggregate(
[
  { $match: { "acceptedNameUsage" : "Tachiramantis douglasi"}},
  { $lookup: {
     "from": "records",
     "foreignField": "taxID",
     "localField": "taxID",
     "as": "alias_records"
      }
  },
  {
      $unwind: "$alias_records"
  },
  {
      $project: {
          
          taxVerifSource: "$alias_records.source",
          country: "$alias_records.country" 
          
          
      }
  }
]
)/*habilitar si se quiere obtener el conteo .itcount()*/