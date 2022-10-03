# erkin
Solução, resolução (gabarito) de todos os problemas (exercícios) do método diferencial, integral de cinética.



<h1><div ;="" align="center">
                 <big style="color: rgb(100, 150, 200);">Ensino de Engenharia Química com Python</big></div>
</h1>
<table align="center">
<tbody>
<tr>
      <td width="300" align="center">
                                    <h1><a href="http://virtualabs.tecnico.ulisboa.pt/virtual/index.html">Virtualabs</a>
                                    </h1></td>
      <td width="100" align="center"><br></td>
      <td width="400" align="center">
              <big style="color: rgb(91, 179, 229);"><h1>Python</h1></big></td>
</tr>
<tr>
      <td align="center">
             <a href="http://custodians.online/portuguese.html">En solidaridad con Library Genesis y Sci-Hub</a>
                                </td>
      <td><br></td>
      <td align="center">SAGE e Python em Engenharia Químicas</td>
</tr>
<tr>
      <td align="center">
            <a href="https://gamafreire.github.io">VIRTUAis LABoratoriaiS no GitHub</a></td>
       <td><br></td>
      <td align="center">Python em português</td></tr>
<tr>
      <td width="10" align="center"><a href="http://visual.icse.us.edu.pl/methodology/sage_in_nutshell.html">SageMath in nutshell</a></td>
      <td align="center"><a>
            <img style="border: 0px solid ;" alt="OpenSource" src="http://virtualabs.ist.utl.pt/osi_keyhole.png" align="middle"></a></td>
        <td align="center">Incrementar o uso de OpenSource no ensino
             <a><img style="border: 0px solid ;" alt="Python" src="http://virtualabs.ist.utl.pt/python-logo_4.png" align="middle"></a>
</td></tr>
</tbody>
</table>
                    
<div class="jp-RenderedHTMLCommon jp-RenderedMarkdown jp-MarkdownOutput " data-mime-type="text/markdown">
        <p></p>
        <hr style="border: 0pt none ; background-color: rgb(83, 79, 154); height: 5px;">
    </div>
                    
## ERKin

Como confirmar a resolução dos, quase todos, problemas de cinética química de disciplinas como engenharia das reacções ou química-fisíca?

Calculam-se, de forma muito simples para o utilizador, k, n e Co ou ko e Ea de dados cinéticos.

Estas pequenas funções permitem assim resolver a maioria dos problemas de determinação de parâmetros a partir de ensaios em descontínuo.

**ERKin** contem duas funções:
- cinetica  - determina Co, k e n
- ERA - determina ko e Ea


#### Limitações

- **cinetica**                

Os dados (método diferencial ou integral) têm de ser duas listas, uma de concentrações ou conversão obtidos num reactor descontínuo (batch) a temperatura constante e a outra, uma lista com os tempos correspondentes.

A lei de velocidade da reacção tem de ser da forma

$$\Re=kC^n$$

que integrada é ajustada, sem linearização, em ordem a k, Co e n.

Se os dados forem de concentração de produto, têm de ser postos em termos de concentração de reagente.

    
- **ERA**

Os dados têm de ser duas listas. Uma com as temperaturas (ºC ou K) e a outra as constantes de velocidade correspondentes.

Trata-se de ajustar, sem linearizar, a **equação de Arrhenius** aos dados $$ k=ko\ e^{-\dfrac{Ea}{RT}}$$


Tem de fornecer pelo menos 2 constantes em função da temperatura em ºC. Se estiverem em K é só ver exemplo abaixo.

                    
<div class="jp-RenderedHTMLCommon jp-RenderedMarkdown jp-MarkdownOutput " data-mime-type="text/markdown">
        <p></p>
        <hr style="border: 0pt none ; background-color: rgb(83, 79, 154); height: 5px;">
    <a href="./cine.html"> Voltar ao início</a></div>
    
### 1 - Instalação

Fazer download do ficheiro <a href="./erkin.py" download=" erkin.py">erkin.py</a> (v2.0, separador decimal com vírgulas em vez de ponto) e movê-lo para o sítio dos seus cálculos.

Claro que precisa ter instalados os pacotes matplotlib e scipy

### 2 - Testes

Se copiar o ficheiro erkin para uma célula em Jupyter e mandar executar deve obter o mesmo resultado que:
