import torch
import torch.nn as nn
import torch.nn.functional as F
import dhg

from dhg.structure.graphs import Graph
from sklearn.metrics.pairwise import cosine_similarity
from typing import List

class HGNNPredictor(nn.Module):
    """The HGNN :sup:`+` model proposed in `HGNN+: General Hypergraph Neural Networks <https://ieeexplore.ieee.org/document/9795251>`_ paper (IEEE T-PAMI 2022).

    Args:
        ``in_channels`` (``int``): :math:`C_{in}` is the number of input channels.
        ``hid_channels`` (``List[int]``): the dims of the HGNNP layer.
        ``out_channels`` (``List[int]``): the dims of the MLP.
        ``use_bn`` (``bool``): If set to ``True``, use batch normalization. Defaults to ``False``.
        ``drop_rate`` (``float``, optional): Dropout ratio. Defaults to ``0.5``.
    """

    def __init__(
        self,
        in_channels: int,
        hgnn_channels: List[int],
        linear_channels: List[int],
        use_bn: bool = None,
        drop_rate: float = 0.5,
    ) -> None:
        super().__init__()
        self.encoder_layers = nn.ModuleList()

        hgnn_channel_list = [in_channels] + hgnn_channels
        for _idx in range(1, len(hgnn_channel_list) - 1):
            self.encoder_layers.append(HGNNP_layer(hgnn_channel_list[_idx-1], hgnn_channel_list[_idx], use_bn=use_bn, drop_rate=drop_rate))

        self.encoder_layers.append(HGNNP_layer(hgnn_channel_list[-2], hgnn_channel_list[-1], use_bn=use_bn, drop_rate=drop_rate,is_last=True))

        self.rec_MLP = MLP_layer(layer_sizes = [hgnn_channel_list[-1]]+linear_channels,use_bn=True)
        self.tf_MLP = MLP_layer(layer_sizes = [hgnn_channel_list[-1]]+linear_channels,use_bn=True)
        self.tg_MLP = MLP_layer(layer_sizes = [hgnn_channel_list[-1]]+linear_channels,use_bn=True)

        # self.linear_decoder = nn.Linear(3*out_channels[-1], 3)

    def encode(self, x: torch.Tensor, hg: "dhg.Hypergraph") -> torch.Tensor:
        """The forward function.
        Args:
            ``x`` (``torch.Tensor``): Input vertex feature matrix. Size :math:`(N, C_{in})`.
            ``hg`` (``dhg.Hypergraph``): The hypergraph structure that contains :math:`N` vertices.
        """

        for layer in self.encoder_layers:
            x = layer(x, hg)
        x = F.elu(x)
        return x

    def decode(self, rec_out: torch.Tensor, tf_out: torch.Tensor, tg_out: torch.Tensor):

        prob1 = torch.cosine_similarity(rec_out,tf_out,dim=1)
        prob2 = torch.cosine_similarity(tf_out,tg_out,dim=1)
        prob = torch.abs(prob1*prob2)
        return prob.view(-1,1)

        # h = torch.cat([rec_out, tf_out, tg_out],dim=1)
        # prob = self.linear_decoder(h)

        # return prob     

    def forward(self, x: torch.Tensor, hg: "dhg.Hypergraph", training_set):

        all_embeddings = self.encode(x, hg)

        
        rec_embd = self.rec_MLP(all_embeddings)
        tf_embd = self.tf_MLP(all_embeddings)
        tg_embd = self.tg_MLP(all_embeddings)

        self.rec_output = rec_embd[training_set[:, 0]]
        self.tf_output = tf_embd[training_set[:, 1]]
        self.tg_output = tg_embd[training_set[:, 2]]   

        """
        rec_embd = all_embeddings[training_set[:, 0]]
        tf_embd = all_embeddings[training_set[:, 1]]
        tg_embd = all_embeddings[training_set[:, 2]]
        self.rec_output = self.rec_MLP(rec_embd)
        self.tf_output = self.tf_MLP(tf_embd)
        self.tg_output = self.tg_MLP(tg_embd)
        """

        return self.decode(self.rec_output, self.tf_output, self.tg_output)
    
    """
    def reset_parameters(self):
        for layer in self.encoder_layers:
            nn.init.xavier_uniform_(layer.weight.data)

        for layer in self.MLP_layers:
            nn.init.xavier_uniform_(layer.weight, gain=1.414)
    """

    def get_embedding(self, x: torch.Tensor, hg: "dhg.Hypergraph", training_set):

        return self.rec_output, self.tf_output, self.tg_output
    

class HGNNP_layer(nn.Module):

    # Modified from: 
    """The HGNN :sup:`+` convolution layer proposed in `HGNN+: General Hypergraph Neural Networks <https://ieeexplore.ieee.org/document/9795251>`_ paper (IEEE T-PAMI 2022).

    Sparse Format:
    
    .. math::

        \left\{
            \begin{aligned}
                m_{\beta}^{t} &=\sum_{\alpha \in \mathcal{N}_{v}(\beta)} M_{v}^{t}\left(x_{\alpha}^{t}\right) \\
                y_{\beta}^{t} &=U_{e}^{t}\left(w_{\beta}, m_{\beta}^{t}\right) \\
                m_{\alpha}^{t+1} &=\sum_{\beta \in \mathcal{N}_{e}(\alpha)} M_{e}^{t}\left(x_{\alpha}^{t}, y_{\beta}^{t}\right) \\
                x_{\alpha}^{t+1} &=U_{v}^{t}\left(x_{\alpha}^{t}, m_{\alpha}^{t+1}\right) \\
            \end{aligned}
        \right.

    Matrix Format:

    .. math::
        \mathbf{X}^{\prime} = \sigma \left( \mathbf{D}_v^{-1} \mathbf{H} \mathbf{W}_e 
        \mathbf{D}_e^{-1} \mathbf{H}^\top \mathbf{X} \mathbf{\Theta} \right).

    Args:
        ``in_channels`` (``int``): :math:`C_{in}` is the number of input channels.
        ``out_channels`` (int): :math:`C_{out}` is the number of output channels.
        ``bias`` (``bool``): If set to ``False``, the layer will not learn the bias parameter. Defaults to ``True``.
        ``use_bn`` (``bool``): If set to ``True``, the layer will use batch normalization. Defaults to ``False``.
        ``drop_rate`` (``float``): If set to a positive number, the layer will use dropout. Defaults to ``0.5``.
        ``is_last`` (``bool``): If set to ``True``, the layer will not apply the final activation and dropout functions. Defaults to ``False``.
    """

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        leaky_alpha: float = 0.1,
        bias: bool = True,
        use_bn: bool = True,
        drop_rate: float = 0.5,
        is_last: bool = True,
    ):
        super().__init__()
        self.is_last = is_last
        self.bn = nn.BatchNorm1d(out_channels) if use_bn else None
        self.act_func = nn.LeakyReLU(negative_slope=leaky_alpha)
        self.drop = nn.Dropout(drop_rate)
        self.theta = nn.Linear(in_channels, out_channels, bias=bias)

    def forward(self, x: torch.Tensor, hg: "dhg.Hypergraph") -> torch.Tensor:
        """The encoder function.

        Args:
            x (``torch.Tensor``): Input vertex feature matrix. Size :math:`(|\mathcal{V}|, C_{in})`.
            hg (``dhg.Hypergraph``): The hypergraph structure that contains :math:`|\mathcal{V}|` vertices.
        """
        x = self.theta(x)
        x = hg.v2v(x, aggr="mean")
        if not self.is_last:
            x = self.act_func(x)
            if self.bn is not None:
                x = self.bn(x)
            x = F.dropout(x,p=0.01)
            x = F.normalize(x,p=2,dim=1)

        return x


class MLP_layer(nn.Module):
    """A common Multi-Layer Perception (MLP) model.

    Args:
        ``channel_list`` (``List[int]``): The list of channels of each layer (input, hiddens, output).
        ``act_name`` (``str``): The name of activation function can be any `activation layer <https://pytorch.org/docs/stable/nn.html#non-linear-activations-weighted-sum-nonlinearity>`_ in Pytorch.
        ``act_kwargs`` (``dict``, optional): The keyword arguments of activation function. Defaults to ``None``.
        ``use_bn`` (``bool``): Whether to use batch normalization.
        ``drop_rate`` (``float``): Dropout ratio. Defaults to ``0.5``.
        ``is_last`` (``bool``): If set to True, the last layer will not use activation, batch normalization, and dropout.
    """

    def __init__(self, layer_sizes, use_bn = True):
        super(MLP_layer, self).__init__()
        self.layers = nn.ModuleList()
        self.norms = nn.ModuleList()
        self.use_bn = use_bn
        for i in range(len(layer_sizes) - 1):
            self.layers.append(nn.Linear(layer_sizes[i], layer_sizes[i+1]))
            # Add a BatchNorm1d layer after each Linear layer except for the output layer
            if i < len(layer_sizes) - 2:  # No BatchNorm after the last linear layer
                self.norms.append(nn.BatchNorm1d(layer_sizes[i+1]))

    def forward(self, x):
        for i, layer in enumerate(self.layers[:-1]):
            x = layer(x)
            if self.use_bn:
                x = self.norms[i](x)  # Apply Batch Normalization
            x = F.elu(x)
        x = self.layers[-1](x)  # No BatchNorm and no ReLU before the final layer

        return x