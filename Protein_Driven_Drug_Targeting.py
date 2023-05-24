import os
import json
import textwrap
import pandas as pd
import time
from bs4 import BeautifulSoup
from elasticsearch import Elasticsearch
from pyecharts.charts import Sunburst ,Tree ,Bar ,Page
from pyecharts import options as opts

# 获取当前目录下txt结尾的文件名
def get_txt():
    for file_name in os.listdir():
        if file_name.endswith('.txt'):
            return file_name

# 获取symbol列表
def get_Symbol(file_name):
    with open(file_name, 'r', encoding='utf-8') as f:
        f = f.readlines()
    Symbol = []
    for i in f:
        i = i.replace('\t', ' ').replace('\n', '').split(' ')
        Symbol.append(i[0])
    return Symbol

# 将组学数据中与每个差异表达蛋白相互作用的蛋白进行汇总，统计每个蛋白出现的次数
# 根据interaction_num的值，筛选出与差异表达蛋白相互作用次数大于等于interaction_num的蛋白
def get_targetNum_dict(symbol_list, interaction_num, PPI_DICT):
    ALL_PPI_PROTEIN = []
    for i in symbol_list:
        if i in PPI_DICT.keys():
            ALL_PPI_PROTEIN.extend(PPI_DICT[i])
    PPI_NUMBER = {}
    TARGET_PPI = {}
    for i in ALL_PPI_PROTEIN:
        if i not in PPI_NUMBER.keys():
            PPI_NUMBER[i] = 1
        else:
            PPI_NUMBER[i] += 1
            if PPI_NUMBER[i] >= interaction_num and i not in TARGET_PPI.keys():
                TARGET_PPI[i] = PPI_NUMBER[i]
            elif PPI_NUMBER[i] >= interaction_num and i in TARGET_PPI.keys():
                TARGET_PPI[i] = PPI_NUMBER[i]
    return TARGET_PPI

# 获取PPI列表
def get_PPI_Symbol_List(symbol_list, interaction_num):
    # 1. 获取PPI的字典
    with open ('data/PPI/PPI.json', 'r') as f:
        PPI_DICT = json.load(f)
    # 2. 筛选出相互作用的蛋白
    TARGET_PPI = get_targetNum_dict(symbol_list, interaction_num, PPI_DICT)
    # 3.获得与差异表达蛋白相互作用的蛋白列表, 还需要将本身这些差异表达蛋白去掉
    TARGET_PPI_LIST = [i for i in TARGET_PPI.keys() if i not in symbol_list]
    return TARGET_PPI_LIST

# 获取靶标分类信息的列表
def classify_targets(Symbol_To_Target, Symbol_list):
    target_have_drug,target_no_drug = [],[]
    target_FDA_approved, target_clinical_trial, target_others = [],[],[]
    
    for symbol in Symbol_list:
        if symbol in Symbol_To_Target.keys():
            target_have_drug.append(symbol)
        else:
            target_no_drug.append(symbol)

    for symbol in target_have_drug:
        target_phage = [*Symbol_To_Target[symbol].values()][0]
        target_name = [*Symbol_To_Target[symbol].keys()][0]
        drug_phase, drug_ap_cl, drug_ap, drug_cl = drug_classify(target_name)
        # 增加一个判断，如果symbol所在阶段确实存在对应的药物，则将其加入到对应的列表中 
        if target_phage == 'Successful target' and drug_ap != []:
            target_FDA_approved.append(symbol)
        elif target_phage == 'Clinical Trial target' and drug_cl != []:
            target_clinical_trial.append(symbol)
        else:
            target_others.append(symbol)
    return target_have_drug, target_no_drug, target_FDA_approved, target_clinical_trial, target_others

# 将html中的数据进行更改后输出
def classify_targets_html(target_have_drug, target_no_drug, target_FDA_approved,
                        target_clinical_trial, target_others, dir_name):
    text_html = open(r'Template/target_pie_template.html',
                    'r', encoding='utf-8').read()
    text_html = text_html.replace(
        'Compound data', str(len(target_have_drug))).replace(
        'No-drug data', str(len(target_no_drug))).replace(
        'FDA Approved data', str(len(target_FDA_approved))).replace(
        'Others data', str(len(target_others))).replace(
        'Clinical data', str(len(target_clinical_trial)))

    soup = BeautifulSoup(text_html, 'html.parser')
    with open('output/' + disease_name + ' reported_number_' + str(reported_number) + '/' + dir_name + '/Targets_pie_chart.html', 'w', encoding='utf-8') as fp:
        fp.write(str(soup))

# 查询既要满足摘要在限定的列表中 又要满足这些摘要中存在HCC这个词组 
# 还要在全部的摘要中输入的symbol和uniprotID  
# hepatocellular carcinoma
def query_target(symbol,Symbol_To_PubMedID,Symbol_To_UniprotID,Symbol_To_Fullname,es,keywords):
    uniprotID = Symbol_To_UniprotID[symbol]
    pubMedId = Symbol_To_PubMedID[symbol]
    # final_list = []
    sql1 = {
        'query': {
            'bool': {
                'must': [
                    {
                        'terms': {
                            'pubMedId': pubMedId
                        }
                    },
                    {
                        "match_phrase": {
                            "abstract": keywords
                            }  
                    },
                    {
                        "match_phrase": {  # abstract中还要存在另一个关键词
                            "abstract": symbol
                            }  
                    },
                ]
            }
        }
    }
    res = es.search(index='abstract22', body=sql1, scroll='5m')
    reported_number_1 = res['hits']['total']['value']
    # if res['hits']['total']['value'] > 0:
    #     for i in res['hits']['hits']:
    #         final_list.append(i['_id'])
    #     final_list = final_list[:5]
    # else:
    #     final_list = []
    es.clear_scroll(scroll_id=res['_scroll_id'])
    #在全部的摘要中检索疾病和symbol是否同时出现
    if reported_number_1 == 0:
        if symbol in Symbol_To_Fullname.keys():
            fullName = Symbol_To_Fullname[symbol]
            sql2 = {
                "query": {
                        "bool": {
                            'must': [
                                {
                                    "match_phrase": {
                                        "abstract": keywords
                                        }  
                                },
                                {
                                    "bool": {
                                        "should": [
                                            {
                                                "match_phrase": {  
                                                    "abstract": fullName
                                                }  
                                            },
                                            {
                                                "match_phrase": {  
                                                    "abstract": symbol
                                                }  
                                            }
                                        ]
                                    }
                                }
                            ]
                        }
                    }
            }
            res = es.search(index='abstract22', body=sql2, scroll='5m')
            reported_number_2 = res['hits']['total']['value']
            # if res['hits']['total']['value'] > 0:
            #     for i in res['hits']['hits']:
            #         final_list.append(i['_id'])
            #     final_list = final_list[:5]
            # else:
            #     final_list = []
            es.clear_scroll(scroll_id=res['_scroll_id'])
        else:
        # 如果没有对应的全称，就将reported_number_2设置为0
        # 这里可以还需要改一改，看看这种情况是否需要使用symbol进行检索
            reported_number_2 = 0
    else:
        reported_number_2 = 0
    # return (reported_number_1 + reported_number_2), final_list
    return reported_number_1 + reported_number_2

# 通过摘要中的关键词进行查询，将靶标分为对于该疾病报道过的靶标和没有报道过的靶标
def report_info(fa, ct, keywords, input_num):
    with open ('data\ID_Transformed\Symbol_To_PubMedID.json', 'r') as f:
        Symbol_To_PubMedID = json.load(f)
    with open ('data\ID_Transformed\Symbol_To_UniprotID.json', 'r') as f:
        Symbol_To_UniprotID = json.load(f)    
    with open ('data\ID_Transformed\Symbol_To_Fullname.json', 'r') as f:
        Symbol_To_Fullname = json.load(f)
    fda_no_review, fda_review ,ct_no_review, ct_review= [],[],[],[]
    for symbol in fa:
        # 存在有的symbol没有对应的uniprotID或者pubMedID 对于这样的symbol进行剔除
        if symbol in Symbol_To_PubMedID.keys() and symbol in Symbol_To_UniprotID.keys():
            query_num = query_target(symbol,Symbol_To_PubMedID,Symbol_To_UniprotID,Symbol_To_Fullname,es,keywords)
            if query_num > input_num:
                fda_review.append(symbol)
            else:
                fda_no_review.append(symbol)
    for symbol in ct:
        if symbol  in Symbol_To_PubMedID.keys() and symbol  in Symbol_To_UniprotID.keys():
            query_num = query_target(symbol,Symbol_To_PubMedID,Symbol_To_UniprotID,Symbol_To_Fullname,es,keywords)
            if query_num > input_num:
                ct_review.append(symbol)
            else:
                ct_no_review.append(symbol)
    return fda_no_review, fda_review ,ct_no_review, ct_review

# 生成靶标信息的Tree图
def all_targets_tree(fda_unre,fda_re,clinical_unre,clinical_re,dir_name):
    fda_all_nmb = len(fda_unre) + len(fda_re)
    cli_all_nmb = len(clinical_unre) + len(clinical_re)
    fda_unre_dict = [{'name': i} for i in fda_unre]
    fda_re_dict = [{'name': i} for i in fda_re]
    clinical_unre_dict = [{'name': i} for i in clinical_unre]
    clinical_re_dict = [{'name': i} for i in clinical_re]
    data = [
    {
        "children": [
            
            {
                "children": [
                    {
                        "children": fda_unre_dict,
                        "name": "Not Reported",
                        'value': len(fda_unre)
                    },
                    {
                        "children":fda_re_dict,
                        "name": "Reported",
                        'value': len(fda_re)
                    }],
                "name": "FDA Approved",
                'value': fda_all_nmb,    
            },
            {
                "children": [
                    {
                        "children":clinical_unre_dict,
                        "name": "Not Reported",
                        'value': len(clinical_unre)
                    },
                    {
                        "children":clinical_re_dict,
                        "name": "Reported",
                        'value': len(clinical_re)
                    }
                ],
                "name": "Clinical",
                'value': cli_all_nmb
            },
        ],
        "name": "Target",
        'value': fda_all_nmb + cli_all_nmb
    }
]
    
    c = (
    Tree(init_opts=opts.InitOpts(width="1200px", height="800px",renderer="svg"))
    .add("", data, collapse_interval=2,
        symbol_size=10,
        symbol="emptyCircle",
        leaves_label_opts=opts.LabelOpts(position="right"),
        itemstyle_opts=opts.ItemStyleOpts(border_width=1, border_color="#48466d"),
        )
    .set_global_opts(title_opts=opts.TitleOpts(title="All Targets"))
    .set_series_opts(label_opts=opts.LabelOpts(
                                            font_size=20,
                                            font_weight='bold',
                                            color="#48466d"
                                            ),
                    )
    
    .render('output/' + disease_name + ' reported_number_' + str(reported_number) + '/' + dir_name + "/Targets_tree.html")
)

# 'hepatocellular carcinoma'
# 查询药物关于疾病的报道信息
def get_drug_report_info(drug_ap, drug_cl, disease, input_num):
    drug_ap_not_report, drug_ap_report, drug_cl_not_report, drug_cl_report = [],[],[],[]
    for drug_name in drug_ap:
        # 这个名字需要进行处理
        # 会出现特殊字符无法处理的情况[Avastin+/-Tarceva]
        drug_name = drug_name.replace(
            '+/-', ' ').replace(
                '/', ' ').replace(
                    '[', '').replace(
                        ']', '').replace(
                            '-', ' ')
        query = {
            'query': {
                'bool': {
                    'must': [
                        {
                            "match": {
                                "abstract": drug_name
                            }
                        },
                        {
                            "match_phrase": {
                                "abstract": disease
                                }  
                        },
                    ]
                }
            }
        }
        res = es.search(index='abstract22', body=query, scroll='5m')
        reported_number = res['hits']['total']['value']
        es.clear_scroll(scroll_id=res['_scroll_id'])
        if reported_number > input_num:
            drug_ap_report.append(drug_name)
        else:
            drug_ap_not_report.append(drug_name)
    for drug_name in drug_cl:
        drug_name = drug_name.replace(
            '+/-', ' ').replace(
                '/', ' ').replace(
                    '[', '').replace(
                        ']', '').replace(
                            '-', ' ')
        query = {
            'query': {
                'bool': {
                    'must': [
                        {
                            "match": {
                                "abstract": drug_name
                            }
                        },
                        {
                            "match_phrase": {
                                "abstract": disease
                                }  
                        },
                    ]
                }
            }
        }
        res = es.search(index='abstract22', body=query, scroll='5m')
        reported_number = res['hits']['total']['value']
        es.clear_scroll(scroll_id=res['_scroll_id'])
        if reported_number > input_num:
            drug_cl_report.append(drug_name)
        else:
            drug_cl_not_report.append(drug_name)
    return drug_ap_not_report, drug_ap_report, drug_cl_not_report, drug_cl_report

# 通过phase将药物进行分类
def drug_classify(target_name):
    with open ('data\Drug\Target_To_Drug.json', 'r') as f:
        Target_To_Drug = json.load(f)
    
    drug_phase = {'Approved': [], 'Clinical_trial': [], 'Others': []}
    drug_ap_cl, drug_ap, drug_cl = [],[],[]
    for drug in Target_To_Drug[target_name]:
        value = [*drug.values()][0]
        key = [*drug.keys()][0]
        if value == 'Approved':
            drug_phase['Approved'].append({
                'name': key
            })
            drug_ap_cl.append(key)
            drug_ap.append(key)
        elif value.startswith('Phase') or value.startswith('Clinical'):
            drug_phase['Clinical_trial'].append({
                'name': key
            })
            drug_ap_cl.append(key)
            drug_cl.append(key)
        else:
            drug_phase['Others'].append({
                'name': key
            })
    return drug_phase, drug_ap_cl, drug_ap, drug_cl

# 从药物列表获取药物的频率
def get_drug_frequency(drug_not_report, drug_report):
    drug_frequency = []
    if drug_not_report!=[]:
        for drug in drug_not_report:
                
            # 这个名字需要进行处理
            # 会出现特殊字符无法处理的情况[Avastin+/-Tarceva]
            drug_name = drug.replace(
                '+/-', ' ').replace(
                    '/', ' ').replace(
                        '[', '').replace(
                            ']', '').replace(
                                '-', ' ')
            query = {
                'query': {
                        "match": {
                            "abstract": drug_name
                            }
                }
            }
            res = es.search(index='abstract22', body=query, scroll='5m')
            reported_number = res['hits']['total']['value']
            drug_frequency.append(reported_number)  
            es.clear_scroll(scroll_id=res['_scroll_id'])
    else:
        for drug in drug_report:
            drug_name = drug.replace(
                '+/-', ' ').replace(
                    '/', ' ').replace(
                        '[', '').replace(
                            ']', '').replace(
                                '-', ' ')
            query = {
                'query': {
                        "match": {
                            "abstract": drug_name
                            }
                }
            }
            res = es.search(index='abstract22', body=query, scroll='5m')
            reported_number = res['hits']['total']['value']
            drug_frequency.append(reported_number)  
            es.clear_scroll(scroll_id=res['_scroll_id'])
    return drug_frequency

# 将字符串处理，过长的字符串换行
def wrap_text(text, max_length=20):
    if len(text) > max_length:
        return textwrap.fill(text, max_length)
    return text

# 将药物转化为tree图需要的数据类型
def drug_treetype_data(drugs):
    drug = []
    for i in drugs:
        name = wrap_text(i)
        drug.append({'name': name})
    return drug

# 生成每个靶标对应药物的树状图以及柱状图
def target_tree_bar(dir_name,symbol,drug_frequency,
                    drug_ap_not_report,drug_ap_report,
                    drug_cl_not_report,drug_cl_report):
    drug_ap_cl = drug_ap_not_report + drug_ap_report + drug_cl_not_report + drug_cl_report
    drug_not_report = drug_ap_not_report + drug_cl_not_report
    drug_report = drug_ap_report + drug_cl_report
    drug_ap_not_report = drug_treetype_data(drug_ap_not_report)
    drug_ap_report = drug_treetype_data(drug_ap_report)
    drug_cl_not_report = drug_treetype_data(drug_cl_not_report)
    drug_cl_report = drug_treetype_data(drug_cl_report)
    tree_data = [
    {
        "children": [
            
            {
                "children": [
                    {
                        "children": drug_ap_not_report,
                        "name": "Not Reported",
                        'value': len(drug_ap_not_report)
                    },
                    {
                        "children":drug_ap_report,
                        "name": "Reported",
                        'value': len(drug_ap_report)
                    }],
                "name": "FDA Approved",
                'value': len(drug_ap_not_report + drug_ap_report),    
            },
            {
                "children": [
                    {
                        "children":drug_cl_not_report,
                        "name": "Not Reported",
                        'value': len(drug_cl_not_report)
                    },
                    {
                        "children":drug_cl_report,
                        "name": "Reported",
                        'value': len(drug_cl_report)
                    }
                ],
                "name": "Clinical",
                'value': len(drug_cl_report + drug_cl_not_report)
            },
        ],
        "name": symbol,
        'value': len(drug_ap_cl)
    }
]

    tree = (
    Tree(init_opts=opts.InitOpts(width="1600px", height="800px",renderer="svg"))
    .add("", tree_data, collapse_interval=2,
        symbol_size=10,
        symbol="emptyCircle",
        leaves_label_opts=opts.LabelOpts(position="right"),
        itemstyle_opts=opts.ItemStyleOpts(border_width=1, border_color="#48466d"),
        edge_fork_position="100%",
        )
    .set_global_opts(title_opts=opts.TitleOpts(title=symbol + ' corresponding drugs'))
    .set_series_opts(label_opts=opts.LabelOpts(
                                            font_size=15,
                                            font_weight='bold',
                                            color="#48466d"
                                            ),
                    )
)
    if drug_not_report == []:
        drug_data = drug_report
    else:
        drug_data = drug_not_report
        
    bar = (
    Bar(init_opts=opts.InitOpts(width="2000px", height="800px",renderer="svg"))
    .add_xaxis(drug_data)
    .add_yaxis("Drug Frequency", drug_frequency, color='#617bdb' )
    .set_global_opts(
        xaxis_opts=opts.AxisOpts(
                                is_show=False,
                                axislabel_opts=opts.LabelOpts(font_size=15,
                                        font_weight='bold',
                                        color="#48466d",
                                        ),
),
        yaxis_opts=opts.AxisOpts(
            
            axislabel_opts=opts.LabelOpts(font_size=15,
                                        font_weight='bold',
                                        color="#48466d"),
            
        ),
        legend_opts=opts.LegendOpts(is_show=False),
        
    )
    .set_series_opts(
        label_opts=opts.LabelOpts(position="right",
                                color="#48466d",
                                font_size=15,
                                font_weight='bold',),
    )
    .reversal_axis()
    )
    
    (
        Page()
        .add(tree,bar)
        ).render(
            'output/' + disease_name + ' reported_number_' + str(reported_number) + '/' + dir_name + '/' + symbol + '/'+ symbol + '.html'
            )

# 制作sunburst图
def get_sunburst(un_relevant_targets_recommend_drug, fa, dir_name):
    data2 = [
        {
            'name': 'FDA approve',
            "itemStyle": {"color": '#fac858'},
            'children': []

        },
            {
            'name': 'Clinical trial',
            "itemStyle": {"color": '#73c0de'},
            'children': []

        },
    ]

    for key, value in un_relevant_targets_recommend_drug.items():
        if key in fa:
            
            children = {
                "name": key,
                'value': 1,
                "itemStyle": {"color":'#fac858' },
                'children': [
                        {'name': value,
                        'value': 1,
                        "itemStyle": {"color":'#fac858' },
                        }
                ]
            }
            data2[0]['children'].append(children)    
        else:
            children = {
                "name": key,
                'value': 1,
                "itemStyle": {"color": '#73c0de'},
                'children': [
                        {'name': value,
                        'value': 1,
                        "itemStyle": {"color": '#73c0de'}
                        }
                ]
            }
            data2[1]['children'].append(children)


    c = (
        Sunburst(init_opts=opts.InitOpts(width="1200px", height="1200px",renderer="svg"))
        .add(
            "",
            data_pair=data2,
            highlight_policy="ancestor",
            sort_="null",
            radius=[0, "95%"],
            # center=["55%", "55%"],
            # 居中
            
            levels=[
                {},
                {
                    "r0": "15%",
                    "r": "35%",
                    "itemStyle": {"borderWidth": 2},
                    "label": {"rotate": "tangential",},
                },
                {"r0": "35%", "r": "60%", "label": {"align": "right"}},
                {
                    "r0": "60%",
                    "r": "62%",
                    "label": {"position": "outside", "padding": 3, "silent": False},
                    "itemStyle": {"borderWidth": 3},
                },
            ],
        )
        .set_global_opts(title_opts=opts.TitleOpts(title="Suggestions for drug targets",
                                                    pos_left='center',
                                                    ))
        .set_series_opts(label_opts=opts.LabelOpts(formatter="{b}",
                                                color='black',
                                                font_weight='bold',
                                                font_size=15,
                                                font_family='Microsoft YaHei',                   
                                                ))
        
        .render('output/' + disease_name + ' reported_number_' + str(reported_number) + '/' + dir_name + "/drug_suggestion.html")
    )

# 生成excel表格
def get_excel(un_relevant_targets_recommend_drug, dir_name):
    # 生成excel表格
    df = pd.DataFrame(un_relevant_targets_recommend_drug.items(), columns=['Target', 'Drug'])
    df.to_excel('output/' + disease_name + ' reported_number_' + str(reported_number) + '/' + dir_name + "/drug_suggestion.html", index=False)

# 生成靶标对应药物的sunburst图和每个靶标对应的药物信息
def get_sunburst_tree_bar(dir_name, fda_no_review, ct_no_review, fa, disease, input):
    target_not_report = fda_no_review + ct_no_review
    un_relevant_targets_recommend_drug = {}
    for symbol in target_not_report:
        target = [*Symbol_To_Target[symbol].keys()][0]
        # print(symbol)
    # 一个靶标对应的药物信息
        drug_phase, drug_ap_cl, drug_ap, drug_cl = drug_classify(target)

        (drug_ap_not_report,
            drug_ap_report, 
            drug_cl_not_report, # 需要提供两个关键词 1.疾病的名字 2.命中的数量
            drug_cl_report) = get_drug_report_info(drug_ap, drug_cl, disease, input)

    # 药物热度频率
        drug_not_report = drug_ap_not_report + drug_cl_not_report
        drug_report = drug_ap_report + drug_cl_report
        drug_frequency = get_drug_frequency(drug_not_report, drug_report)
        if drug_frequency != []:
            os.makedirs('output/' + disease_name + ' reported_number_' + str(reported_number) + '/' + dir_name + '/'+ symbol, exist_ok=True)
            target_tree_bar(dir_name, symbol,drug_frequency,
                                drug_ap_not_report,drug_ap_report,
                                drug_cl_not_report,drug_cl_report)

            number_index = drug_frequency.index(max(drug_frequency)) 
            if drug_not_report == []:
                suggest_drug = drug_report[number_index]
            else:
                suggest_drug = drug_not_report[number_index]
            un_relevant_targets_recommend_drug[symbol] = suggest_drug
            
    # 输出为excel文件       
    df = pd.DataFrame(un_relevant_targets_recommend_drug.items(), columns=['Target', 'Recommend Drug'])
    df.to_excel('output/' + disease_name + ' reported_number_' + str(reported_number) + '/' + dir_name + "/drug_suggestion.xlsx", index=False)
    get_sunburst(un_relevant_targets_recommend_drug, fa, dir_name)

# 对推荐靶标数量进行控制
def sort_targets(no_review, target_max_number):
    # 将靶标和对应文献数量做成列表
    sort_list = []
    for target in no_review:
        res = es.search(index="abstract22", body={"query": {"match": {"abstract": target}}},scroll='5m')
        target_hot = res['hits']['total']['value']
        es.clear_scroll(scroll_id=res['_scroll_id'])
        sort_list.append([target, target_hot])

    # 使用sorted函数对列表进行排序
    sort_list = sorted(sort_list, key=lambda x: x[1], reverse=True)
    sort_list = [x[0] for x in sort_list]

    # 判断推荐的靶标数量是否大于靶标推荐最大值
    if len(sort_list) > target_max_number:
        sort_list = sort_list[:target_max_number]
    else:
        pass
    return sort_list

# 生成新的靶标列表
def new_targets_list(list,sort_list):
    new_list = [x for x in list if x in sort_list]
    return new_list

# 创建输出的文件夹,基于疾病关键词和reported number
# 获取config文件中的参数
with open('config.json', 'r') as f:
    config = json.load(f)
    
print('Welcome to use the new tool OTTM')
print("-------------------------------------------------------------------")
print('Please make sure to correctly fill in the configuration file')

disease_name = config['disease_name']
reported_number = config['reported_number']
os.makedirs('output/' + disease_name + ' reported_number_' + str(reported_number), exist_ok=True)
os.makedirs('output/' + disease_name + ' reported_number_' + str(reported_number) + '/Target', exist_ok=True)
os.makedirs('output/' + disease_name + ' reported_number_' + str(reported_number) + '/PPI_Target', exist_ok=True)
# 获得组学数据文件的名称
file_name = get_txt()

# 1.获得文件中的symbol
Symbol_list = get_Symbol(file_name)

print('The number of proteins in the file is: ', len(Symbol_list))

# 2.获取symbol对应的PPI蛋白 在这里只要跟差异表达蛋白拥有相互作用就进行保留
# interaction_num = 0 已经存在相互作用的蛋白 最少应该和几个差异表达蛋白列表中的蛋白相互作用 
# 可以设置为参数
interaction_num = config['interaction_num']
Symbol_PPI_list = get_PPI_Symbol_List(Symbol_list, interaction_num)

print('The number of PPI proteins list is: ', len(Symbol_PPI_list))

print('Please wait for a while, the program is running...')

# 显示实时运行时间
t0=time.time()
print('Program start time:',time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))

# 3.获取symbol对应靶标信息
with open ('data/Drug/Symbol_To_Target.json', 'r') as f:
    Symbol_To_Target = json.load(f)

with open ('data/ID_Transformed/Symbol_To_Fullname.json', 'r') as f:
    Symbol_To_Fullname = json.load(f)

# 3.1 PPI对应的靶标分类
p_h_dr,p_no_dr,p_fa,p_ct,p_ot = classify_targets(Symbol_To_Target, Symbol_PPI_list)
# 3.2 symbol对应的靶标分类
h_dr,no_dr,fa,ct,ot = classify_targets(Symbol_To_Target, Symbol_list)

# 3.3 将靶标分类信息输出到html中制作两张图分别保存到Target/PPI_Target文件夹中
classify_targets_html(h_dr,no_dr,fa,ct,ot,'Target')
classify_targets_html(p_h_dr,p_no_dr,p_fa,p_ct,p_ot,'PPI_Target')

# 进行查询
es = Elasticsearch(timeout=30, max_retries=10, retry_on_timeout=True)

# es = Elasticsearch(['https://localhost:9200'], 
#                    http_auth=('elastic', 'zjFz0KjwlZF-2twco1jg'),
#                    verify_certs=False)

# symbol对应的靶标
fda_no_review, fda_review ,ct_no_review, ct_review = report_info(fa,ct,disease_name, reported_number)
p_fda_no_review, p_fda_review ,p_ct_no_review, p_ct_review = report_info(p_fa,p_ct,disease_name, reported_number)

# 生成靶标信息的Tree图
all_targets_tree(fda_no_review, fda_review ,ct_no_review, ct_review,'Target')
all_targets_tree(p_fda_no_review, p_fda_review ,p_ct_no_review, p_ct_review,'PPI_Target')

t1=time.time()
print('Program end time:',time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
print("Running Time：%.6fs"%(t1-t0))



print('The number of recommend targets is: ', len(fda_no_review + ct_no_review))
print('The number of recommend PPI targets is: ', len(p_fda_no_review + p_ct_no_review))

# 生成靶标列表
target_max_number = config['target_max_number'] # 靶标推荐最大值,根据文献数量进行排序

no_review = fda_no_review + ct_no_review
p_no_review = p_fda_no_review + p_ct_no_review

fda_no_review = new_targets_list(fda_no_review,sort_targets(no_review, target_max_number))
ct_no_review = new_targets_list(ct_no_review,sort_targets(no_review, target_max_number))
p_fda_no_review = new_targets_list(p_fda_no_review,sort_targets(p_no_review, target_max_number))
p_ct_no_review = new_targets_list(p_ct_no_review,sort_targets(p_no_review ,target_max_number))

# 获得全部靶标药物推荐的旭日图，以及每个靶标对应的药物信息和药物热度
get_sunburst_tree_bar('Target', fda_no_review, ct_no_review, fa, disease_name, reported_number)
get_sunburst_tree_bar('PPI_Target', p_fda_no_review, p_ct_no_review, p_fa, disease_name, reported_number)

print('The program is finished!')
print('-------------------------------------------------------------------')
input('Please enter any key to exit')

#TODO: 将这个脚本转换成exe文件


# query = {
#     'query': {
#             "match_phrase": {
#                 "abstract": 'hepatocellular carcinomas'
#                 }
#     }
# }
# res = es.search(index='abstract22', body=query, scroll='5m')
# reported_number = res['hits']['total']['value']
# print(reported_number)
# #%%
# def list_info(fa, keywords, input_num, es):
#     list_re = {}
#     list_not_re ={}
#     def query_target_2(symbol, es, keywords):
#         query = {
#             'query': {
#                 'bool': {
#                     'must': [
#                         {
#                             "match": {
#                                 "abstract": symbol
#                             }
#                         },
#                         {
#                             "match_phrase": {
#                                 "abstract": keywords
#                                 }  
#                         },
#                     ]
#                 }
#             }
#         }
#         res = es.search(index='abstract22', body=query, scroll='5m')
#         final_list = []
#         if res['hits']['total']['value'] > 0:
#             for i in res['hits']['hits']:
#                 final_list.append(i['_id'])
#             final_list = final_list[:5]
#         reported_number = res['hits']['total']['value']
#         es.clear_scroll(scroll_id=res['_scroll_id'])
#         return reported_number, final_list
#     with open ('data\ID_Transformed\Symbol_To_PubMedID.json', 'r') as f:
#         Symbol_To_PubMedID = json.load(f)
#     with open ('data\ID_Transformed\Symbol_To_UniprotID.json', 'r') as f:
#         Symbol_To_UniprotID = json.load(f)    
#     with open ('data\ID_Transformed\Symbol_To_Fullname.json', 'r') as f:
#         Symbol_To_Fullname = json.load(f)
#     for symbol in fa:
#         # 存在有的symbol没有对应的uniprotID或者pubMedID 对于这样的symbol进行剔除
#         if symbol in Symbol_To_PubMedID.keys() and symbol in Symbol_To_UniprotID.keys():
#             query_num, Id_list = query_target(symbol,Symbol_To_PubMedID,Symbol_To_UniprotID,Symbol_To_Fullname,es,keywords)
#             if query_num > input_num:
#                 list_re[symbol] = {'reported number': query_num, 'PubMedID': Id_list}
#             else:
                
#                 list_not_re[symbol] = {'reported number': query_num, 'PubMedID': Id_list}
#         else:
#             query_num, Id_list = query_target_2(symbol, es, keywords)
#             if query_num > input_num:
#                 list_re[symbol] = {'reported number': query_num, 'PubMedID': Id_list}
#             else:
#                 print(query_num)
#                 list_not_re[symbol] = {'reported number': query_num, 'PubMedID': Id_list}
#     return list_re, list_not_re

# #%%
# es = Elasticsearch(timeout=30, max_retries=10, retry_on_timeout=True)
# list_re, list_not_re = list_info(h_dr, disease_name, reported_number, es)
# ## no_dr , h_dr

# #%%
# def classify_targets(Symbol_To_Target, Symbol_list):
#     target_have_drug,target_no_drug = [],[]
#     target_FDA_approved, target_clinical_trial, target_others = [],[],[]
#     phase_dict = {}

#     for symbol in Symbol_list:
#         if symbol not in phase_dict.keys():
#             phase_dict[symbol] = {}
#         target_phage = [*Symbol_To_Target[symbol].values()][0]
#         target_name = [*Symbol_To_Target[symbol].keys()][0]
#         drug_phase, drug_ap_cl, drug_ap, drug_cl = drug_classify(target_name)
#         # 增加一个判断，如果symbol所在阶段确实存在对应的药物，则将其加入到对应的列表中 
#         if target_phage == 'Successful target' and drug_ap != []:
#             phase_dict[symbol] = 'FDA Approved Drugs'
#         elif target_phage == 'Clinical Trial target' and drug_cl != []:
#             phase_dict[symbol] = 'Clinical Trial Drugs'
#         else:
#             phase_dict[symbol] = 'Others'
#     return phase_dict
# #%%
# phase_dict = classify_targets(Symbol_To_Target, list_not_re)
# #%%
# df = pd.DataFrame([phase_dict]).T
# df.to_excel('eported.xlsx')