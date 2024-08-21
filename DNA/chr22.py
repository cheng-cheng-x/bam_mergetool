import os
from pyfaidx import Fasta 
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader
from transformers import AdamW

# 定义自定义BERT模型
class CustomBERTForSequenceClassification(nn.Module):
    def __init__(self, num_labels, vocab_size, embed_size):
        super(CustomBERTForSequenceClassification, self).__init__()
        self.embedding = nn.Embedding(vocab_size, embed_size, padding_idx=0)
        self.encoder_layer = nn.TransformerEncoderLayer(d_model=embed_size, nhead=8)
        self.transformer_encoder = nn.TransformerEncoder(self.encoder_layer, num_layers=6)
        self.classifier = nn.Linear(embed_size, num_labels)

    def forward(self, input_ids):
        embeds = self.embedding(input_ids)
        transformer_output = self.transformer_encoder(embeds)
        cls_output = transformer_output[:, 0, :]
        logits = self.classifier(cls_output)
        return logits

# 下载和提取chr22序列
fasta_file = "/data/haocheng/data/DNA/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
print("正在加载FASTA文件...")
fa = Fasta(fasta_file)

print("正在提取chr22序列...")
chr22_sequence = str(fa["22"])  # 提取chr22序列
print("提取完成。序列长度:", len(chr22_sequence))

# 将DNA序列转换为tokens
def dna_to_tokens(sequence):
    base_to_token = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    return [base_to_token[base] for base in sequence]

print("正在将DNA序列转换为tokens...")
tokens = dna_to_tokens(chr22_sequence)  # 转换为tokens
print("转换完成。Tokens示例:", tokens[:100])

# 分割序列
def split_sequence(sequence, chunk_size):
    return [sequence[i:i + chunk_size] for i in range(0, len(sequence), chunk_size)]

chunk_size = 256  # 每个片段包含256个碱基对
print(f"正在将序列分割为大小为 {chunk_size} 的片段...")
chunks = split_sequence(chr22_sequence, chunk_size)  # 分割序列
print("片段创建完成。总片段数:", len(chunks))

# 生成标签
base_to_token = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
labels = [[base_to_token[base] for base in chunk] for chunk in chunks]  # 将每个片段的标签转换为数字
print("标签生成完成。示例标签:", labels[:5])  # 输出前5个标签示例

# 创建数据集
class DNADataset(Dataset):
    def __init__(self, encodings, labels):
        self.encodings = encodings
        self.labels = labels

    def __getitem__(self, idx):
        item = {key: val[idx] for key, val in self.encodings.items()}
        item['labels'] = torch.tensor(self.labels[idx], dtype=torch.long)  # 确保标签为LongTensor
        return item

    def __len__(self):
        return len(self.labels)

    @staticmethod
    def collate_fn(batch):
        input_ids = [item['input_ids'] for item in batch]
        labels = [item['labels'] for item in batch]

        # 填充 input_ids 以确保它们具有相同的长度
        input_ids_padded = nn.utils.rnn.pad_sequence(input_ids, batch_first=True)
        
        # 找到最大长度进行填充
        max_length = input_ids_padded.size(1)
        labels_padded = [torch.cat([label, torch.tensor([4] * (max_length - len(label)), dtype=torch.long)]) for label in labels]  # 填充标签
        labels_padded = torch.stack(labels_padded)  # 将填充后的标签堆叠成一个tensor

        return {'input_ids': input_ids_padded, 'labels': labels_padded}

# 将DNA片段的tokens转换为Tensor并创建数据集
print("正在创建数据集...")
encodings = {'input_ids': [torch.tensor(dna_to_tokens(chunk)) for chunk in chunks]}  # 创建输入ids
dataset = DNADataset(encodings, labels)  # 创建数据集
print("数据集创建完成。总样本数:", len(dataset))

# 创建自定义模型
num_labels = 5  # 设置标签数量
vocab_size = 5  # 4种碱基 + 1个N
embed_size = 128  # 嵌入维度
model = CustomBERTForSequenceClassification(num_labels=num_labels, vocab_size=vocab_size, embed_size=embed_size)  # 使用自定义模型
print("模型初始化完成。")

# 创建DataLoader
batch_size = 8  # 每批次的大小
dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True, collate_fn=DNADataset.collate_fn)  # 创建DataLoader
print("DataLoader创建完成。")

# 设置优化器
optimizer = AdamW(model.parameters(), lr=5e-5)  # 使用AdamW优化器


# 训练模型
num_epochs = 3  # 训练轮数
print("开始训练...")
model.train()  # 设置模型为训练模式
for epoch in range(num_epochs):
    print(f"正在训练第 {epoch + 1}/{num_epochs} 轮...")
    epoch_loss = 0  # 初始化每个epoch的损失
    total_correct = 0  # 初始化正确预测数量
    total_tokens = 0  # 初始化总token数量
    
    for step, batch in enumerate(dataloader):
        # 从数据集中获取输入和标签
        inputs = {key: val for key, val in batch.items()}
        input_ids = inputs.pop("input_ids")
        labels = inputs.pop("labels")

        optimizer.zero_grad()  # 清空梯度
        outputs = model(input_ids)  # 前向传播

        # 计算有效标签的损失，忽略 N 标签
        valid_indices = labels != 4  # 忽略标签为N（即4）

        # 确保 valid_indices 只筛选有效的输出
        if valid_indices.sum() > 0:  # 如果存在有效标签
            # 获取当前批次的有效输出和标签
            valid_outputs = outputs[valid_indices.nonzero(as_tuple=True)[0]]
            valid_labels = labels[valid_indices].long()  # 确保valid_labels是LongTensor

            loss = nn.CrossEntropyLoss()(valid_outputs, valid_labels)  # 计算有效标签的损失
            loss.backward()  # 反向传播
            optimizer.step()  # 更新参数

            # 计算预测的遮住token的准确率
            preds = torch.argmax(valid_outputs, dim=1)
            total_correct += (preds == valid_labels).sum().item()  # 累加正确预测数量
            total_tokens += valid_labels.size(0)  # 累加总token数量

            epoch_loss += loss.item()  # 累加每个batch的损失

            # 打印当前批次的损失、准确率和perplexity
            perplexity = torch.exp(torch.tensor(epoch_loss / total_tokens))  # 计算当前的perplexity
            accuracy = total_correct / total_tokens  # 计算当前的准确率
            print(f"进度: {step}/{len(dataloader)} - 当前批次损失: {loss.item():.4f}，准确率: {accuracy:.4f}，PPL: {perplexity:.4f}")

    # 计算每个epoch的perplexity (PPL) 和准确率
    perplexity = torch.exp(torch.tensor(epoch_loss / total_tokens))  # 使用每个epoch的平均损失计算perplexity
    accuracy = total_correct / total_tokens  # 计算准确率

    print(f"第 {epoch + 1} 轮训练完成。损失: {epoch_loss:.4f}，准确率: {accuracy:.4f}，PPL: {perplexity:.4f}")

print("训练完成。")
